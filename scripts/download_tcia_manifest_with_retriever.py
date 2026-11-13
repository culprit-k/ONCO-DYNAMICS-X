"""Download an explicitly reviewed TCIA SeriesInstanceUID manifest.

The script installs the latest official TCIA Data Retriever Linux x86_64 CLI
release when no executable is supplied, verifies the release asset SHA-256,
and invokes the retriever in conservative/resume mode. It does not accept the
TCIA data-use policy unless the caller explicitly supplies --accept-data-policy.

This is intended for Linux x86_64 Google Colab runtimes. Keep manifests,
logs, DICOMs, and all patient-level files in private storage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

RELEASES_API = "https://api.github.com/repos/TCIA/data-retriever/releases/latest"
USER_AGENT = "ONCO-DYNAMICS-X/0.1 (TCIA public imaging download)"


def _request_json(url: str, timeout: int = 60) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/vnd.github+json"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read().decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Unexpected GitHub release response")
    return payload


def _download_file(url: str, destination: Path, timeout: int = 120) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    with urllib.request.urlopen(request, timeout=timeout) as response, destination.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
            digest.update(chunk)
    return digest.hexdigest()


def _safe_asset_name(release: dict[str, Any]) -> tuple[str, str, str]:
    assets = release.get("assets")
    if not isinstance(assets, list):
        raise ValueError("Latest TCIA Data Retriever release has no assets list")
    candidates = [
        asset for asset in assets
        if isinstance(asset, dict)
        and asset.get("name") == "TCIA_Data_Retriever_linux_x86_64.zip"
        and asset.get("browser_download_url")
    ]
    if len(candidates) != 1:
        raise ValueError(
            "Could not identify exactly one official Linux x86_64 release asset; "
            f"found {len(candidates)}"
        )
    asset = candidates[0]
    digest = str(asset.get("digest") or "")
    if not digest.startswith("sha256:"):
        raise ValueError("Release asset does not publish a SHA-256 digest; refusing unverified binary")
    return (
        str(asset["name"]),
        str(asset["browser_download_url"]),
        digest.split(":", 1)[1].lower(),
    )


def install_latest_retriever(cache_dir: Path) -> tuple[Path, str]:
    """Download and verify the current official Linux x86_64 CLI release."""
    if sys.platform != "linux":
        raise RuntimeError("Automatic installation currently supports Linux/Google Colab only")
    machine = os.uname().machine.lower()
    if machine not in {"x86_64", "amd64"}:
        raise RuntimeError(f"Expected Linux x86_64 runtime; detected {machine}")

    release = _request_json(RELEASES_API)
    version = str(release.get("tag_name") or release.get("name") or "unknown")
    asset_name, asset_url, expected_sha256 = _safe_asset_name(release)
    cache_dir.mkdir(parents=True, exist_ok=True)
    executable = cache_dir / "TCIA_Data_Retriever"
    version_file = cache_dir / "retriever_version.txt"

    if executable.is_file() and os.access(executable, os.X_OK) and version_file.exists():
        cached_version = version_file.read_text(encoding="utf-8").strip()
        if cached_version == version:
            return executable, version

    archive = cache_dir / asset_name
    actual_sha256 = _download_file(asset_url, archive)
    if actual_sha256 != expected_sha256:
        archive.unlink(missing_ok=True)
        raise RuntimeError(
            "TCIA Data Retriever archive SHA-256 verification failed; "
            f"expected {expected_sha256}, got {actual_sha256}"
        )

    with zipfile.ZipFile(archive) as zf:
        candidates = [
            info for info in zf.infolist()
            if not info.is_dir() and Path(info.filename).name == "TCIA_Data_Retriever"
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                "Expected one TCIA_Data_Retriever executable in the official archive; "
                f"found {len(candidates)}"
            )
        info = candidates[0]
        with zf.open(info) as src, executable.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    executable.chmod(0o755)
    version_file.write_text(version + "\n", encoding="utf-8")
    archive.unlink(missing_ok=True)
    return executable, version


def validate_uid_manifest(path: Path) -> list[str]:
    """Validate a one-column SeriesInstanceUID CSV and return unique UIDs."""
    import csv

    if not path.is_file():
        raise FileNotFoundError(f"Manifest does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows or rows[0] != ["SeriesInstanceUID"]:
        raise ValueError("Manifest header must be exactly: SeriesInstanceUID")
    values = [row[0].strip() for row in rows[1:] if row and row[0].strip()]
    if not values:
        raise ValueError("Manifest contains no SeriesInstanceUID values")
    if len(values) != len(set(values)):
        raise ValueError("Manifest contains duplicate SeriesInstanceUID values")
    invalid = [value for value in values if not value.replace(".", "").isdigit()]
    if invalid:
        raise ValueError(f"Manifest contains malformed SeriesInstanceUID values: {invalid[:3]}")
    return values


def run_download(
    manifest: Path,
    output_dir: Path,
    *,
    executable: Path | None = None,
    cache_dir: Path | None = None,
    accept_data_policy: bool = False,
    dry_run: bool = False,
) -> int:
    uids = validate_uid_manifest(manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir or (Path.home() / ".cache" / "onco-dynamics-x" / "tcia-retriever")

    version = "user-supplied"
    if executable is None:
        executable, version = install_latest_retriever(cache)
    else:
        executable = executable.expanduser().resolve()
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise FileNotFoundError(f"Retriever executable is missing or not executable: {executable}")

    command = [
        str(executable),
        "--cli",
        "--input", str(manifest.resolve()),
        "--output", str(output_dir.resolve()),
        "--server-friendly",
        "--skip-existing",
        "--save-log",
    ]
    # This flag represents a user's explicit data-policy acceptance. It is never
    # added automatically by this script.
    if accept_data_policy:
        command.append("--accept-data-policy")

    print(json.dumps({
        "retriever_version": version,
        "manifest": str(manifest.resolve()),
        "series_uid_count": len(uids),
        "output_directory": str(output_dir.resolve()),
        "resume_enabled": True,
        "server_friendly": True,
        "data_policy_accepted_by_caller": bool(accept_data_policy),
        "command": command,
    }, indent=2))
    if dry_run:
        return 0

    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        print(
            f"TCIA Data Retriever exited with status {completed.returncode}. "
            "You can resume with the same command because --skip-existing is enabled.",
            file=sys.stderr,
        )
    return int(completed.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--retriever", type=Path, help="Optional pre-installed executable")
    parser.add_argument(
        "--accept-data-policy",
        action="store_true",
        help="Explicitly pass the TCIA policy-acceptance flag to the retriever; use only after reading and accepting the policy yourself",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    raise SystemExit(run_download(
        args.manifest,
        args.output,
        executable=args.retriever,
        cache_dir=args.cache_dir,
        accept_data_policy=args.accept_data_policy,
        dry_run=args.dry_run,
    ))


if __name__ == "__main__":
    main()
