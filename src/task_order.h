/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_TASK_ORDER_H
#define SWIFT_TASK_ORDER_H

#include "../config.h"

/* Local includes */
#include "scheduler.h"

#ifdef TASK_ORDER_NONE
#include "task_order/none/task_order.h"
#elif TASK_ORDER_GEAR
#include "task_order/GEAR/task_order.h"
#elif TASK_ORDER_EAGLE
#include "task_order/EAGLE/task_order.h"
#else
#error undefined task order
#endif

#endif /* SWIFT_TASK_ORDER_H */
