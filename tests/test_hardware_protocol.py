import pytest

from onco_dynamics_x.hardware.serial_protocol import (
    ControlMessage,
    apply_controls,
    parse_control_message,
)
from onco_dynamics_x.parameters import ModelParameters


def test_parse_valid_control_message():
    message = parse_control_message(
        '{"event":"controls","tumour_growth":0.4,"immune_activity":0.8,"run":true}'
    )
    assert message == ControlMessage(0.4, 0.8, run=True, reset=False)


def test_reject_unknown_control_message_keys():
    with pytest.raises(ValueError, match="Unknown"):
        parse_control_message(
            '{"event":"controls","tumour_growth":0.4,"immune_activity":0.8,"dose":4}'
        )


def test_reject_control_values_outside_normalized_range():
    with pytest.raises(ValueError, match="within"):
        parse_control_message(
            '{"event":"controls","tumour_growth":1.4,"immune_activity":0.8}'
        )


def test_reject_non_boolean_command_flags():
    with pytest.raises(ValueError, match="boolean"):
        parse_control_message(
            '{"event":"controls","tumour_growth":0.4,"immune_activity":0.8,"run":1}'
        )


def test_controls_map_to_bounded_model_parameters():
    message = ControlMessage(tumour_growth=0.5, immune_activity=0.5)
    parameters = apply_controls(ModelParameters(), message)
    assert parameters.cancer_growth == pytest.approx(0.65)
    assert parameters.immune_killing == pytest.approx(0.425)
