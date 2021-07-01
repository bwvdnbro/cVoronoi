/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_HYDRO_SPACE_H
#define SWIFT_HYDRO_SPACE_H

/**
 * @brief Extra space information that is needed for some hydro schemes.
 */
struct hydro_space {
  /*! Anchor of the simulation space. */
  double anchor[3];

  /*! Side lengths of the simulation space. */
  double side[3];
};

inline static void hydro_space_init(struct hydro_space *hs, const double *dim) {
  hs->anchor[0] = 0.;
  hs->anchor[1] = 0.;
  hs->anchor[2] = 0.;
  hs->side[0] = dim[0];
  hs->side[1] = dim[1];
  hs->side[2] = dim[2];
}

#endif /* SWIFT_HYDRO_SPACE_H */
