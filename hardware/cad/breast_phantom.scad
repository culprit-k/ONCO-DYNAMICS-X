// ONCO-DYNAMICS X educational breast phantom
// Units: millimetres. Simplified teaching prop, not anatomical imaging.
part = "shell"; // [shell, small_insert, medium_insert, large_insert, stand_ring]
$fn = 96;
shell_width = 76;
shell_depth = 92;
shell_height = 72;
shell_center_z = 74;
wall_thickness = 3;
opening_height = 13;

module outer_ellipsoid() {
    translate([0, 0, shell_center_z])
        scale([shell_width, shell_depth, shell_height])
            sphere(r = 1);
}
module inner_ellipsoid() {
    translate([0, 0, shell_center_z])
        scale([shell_width - wall_thickness, shell_depth - wall_thickness,
               shell_height - wall_thickness])
            sphere(r = 1);
}
module phantom_shell() {
    difference() {
        difference() {
            outer_ellipsoid();
            inner_ellipsoid();
        }
        translate([-110, -110, -100])
            cube([220, 220, 100 + opening_height]);
    }
}
module lesion_insert(radius = 8) {
    union() {
        sphere(r = radius);
        translate([0, 0, -radius - 7])
            cylinder(h = 7, r = 2.2);
    }
}
module stand_ring() {
    difference() {
        cylinder(h = 7, r = shell_width + 5);
        translate([0, 0, -0.5])
            cylinder(h = 8, r = shell_width - 2);
    }
}
if (part == "shell") {
    color([0.9, 0.72, 0.64]) phantom_shell();
} else if (part == "small_insert") {
    color([0.75, 0.12, 0.18]) lesion_insert(5);
} else if (part == "medium_insert") {
    color([0.75, 0.12, 0.18]) lesion_insert(9);
} else if (part == "large_insert") {
    color([0.75, 0.12, 0.18]) lesion_insert(14);
} else if (part == "stand_ring") {
    color([0.25, 0.28, 0.32]) stand_ring();
}
