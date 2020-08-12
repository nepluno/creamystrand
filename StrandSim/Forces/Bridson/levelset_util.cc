/**
 * \copyright 2008 Robert Bridson
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "levelset_util.hh"

namespace strandsim {
namespace bridson {
// Given two signed distance values (line endpoints), determine what fraction of
// a connecting segment is "inside"
double fraction_inside(double phi_left, double phi_right) {
  if (phi_left < 0 && phi_right < 0) return 1;
  if (phi_left < 0 && phi_right >= 0) return phi_left / (phi_left - phi_right);
  if (phi_left >= 0 && phi_right < 0)
    return phi_right / (phi_right - phi_left);
  else
    return 0;
}

static void cycle_array(double* arr, int size) {
  double t = arr[0];
  for (int i = 0; i < size - 1; ++i) arr[i] = arr[i + 1];
  arr[size - 1] = t;
}

// Given four signed distance values (square corners), determine what fraction
// of the square is "inside"
double fraction_inside(double phi_bl, double phi_br, double phi_tl,
                       double phi_tr) {
  int inside_count = (phi_bl < 0 ? 1 : 0) + (phi_tl < 0 ? 1 : 0) +
                     (phi_br < 0 ? 1 : 0) + (phi_tr < 0 ? 1 : 0);
  double list[] = {phi_bl, phi_br, phi_tr, phi_tl};

  if (inside_count == 4)
    return 1;
  else if (inside_count == 3) {
    // rotate until the positive value is in the first position
    while (list[0] < 0) {
      cycle_array(list, 4);
    }

    // Work out the area of the exterior triangle
    double side0 = 1 - fraction_inside(list[0], list[3]);
    double side1 = 1 - fraction_inside(list[0], list[1]);
    return 1 - 0.5 * side0 * side1;
  } else if (inside_count == 2) {
    // rotate until a negative value is in the first position, and the next
    // negative is in either slot 1 or 2.
    while (list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
      cycle_array(list, 4);
    }

    if (list[1] < 0) {  // the matching signs are adjacent
      double side_left = fraction_inside(list[0], list[3]);
      double side_right = fraction_inside(list[1], list[2]);
      return 0.5 * (side_left + side_right);
    } else {  // matching signs are diagonally opposite
      // determine the centre point's sign to disambiguate this case
      double middle_point = 0.25f * (list[0] + list[1] + list[2] + list[3]);
      if (middle_point < 0) {
        double area = 0;

        // first triangle (top left)
        double side1 = 1 - fraction_inside(list[0], list[3]);
        double side3 = 1 - fraction_inside(list[2], list[3]);

        area += 0.5 * side1 * side3;

        // second triangle (top right)
        double side2 = 1 - fraction_inside(list[2], list[1]);
        double side0 = 1 - fraction_inside(list[0], list[1]);
        area += 0.5 * side0 * side2;

        return 1 - area;
      } else {
        double area = 0;

        // first triangle (bottom left)
        double side0 = fraction_inside(list[0], list[1]);
        double side1 = fraction_inside(list[0], list[3]);
        area += 0.5 * side0 * side1;

        // second triangle (top right)
        double side2 = fraction_inside(list[2], list[1]);
        double side3 = fraction_inside(list[2], list[3]);
        area += 0.5 * side2 * side3;
        return area;
      }
    }
  } else if (inside_count == 1) {
    // rotate until the negative value is in the first position
    while (list[0] >= 0) {
      cycle_array(list, 4);
    }

    // Work out the area of the interior triangle, and subtract from 1.
    double side0 = fraction_inside(list[0], list[3]);
    double side1 = fraction_inside(list[0], list[1]);
    return 0.5 * side0 * side1;
  } else
    return 0;
}
}  // namespace bridson
}  // namespace strandsim
