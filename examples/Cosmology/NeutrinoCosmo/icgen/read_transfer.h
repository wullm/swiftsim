/**
 * @read_transfer.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Read CLASS format transfer function files.
 */

#ifndef READ_TRANSFER_H
#define READ_TRANSFER_H

#include <math.h>
#include <fstream>
#include <vector>

void read_transfer(std::vector<double>&, std::vector<double>&,
                   std::vector<double>&);

#endif
