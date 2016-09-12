/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "LikwidPerformanceMonitoringModule.h"

#ifdef LIKWID_AVAILABLE

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <likwid.h>
#include <sstream>
#include <utility>
#include <vector>

// TODO(guera): Remove
#include <chrono>
#include <thread>

#include "../LikwidProfiler.h"

// TODO(guera): Remove once likwid 4.1 is available on SuperMUC. At the moment
// this section is somewhat Haswell-EP specific.

namespace {
constexpr const int kNumberOfGroups = 16;

constexpr const char* groups[kNumberOfGroups] = {
    "BRANCH",  "CLOCK",    "DATA",      "ENERGY",   "ICACHE", "L2",
    "L2CACHE", "L3",       "L3CACHE",   "MEM",      "NUMA",   "QPI",
    "SBOX",    "TLB_DATA", "TLB_INSTR", "FLOPS_AVX"};
constexpr const char* eventsets[kNumberOfGroups] = {
    /* BRANCH */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "BR_INST_RETIRED_ALL_BRANCHES:PMC0,"
    "BR_MISP_RETIRED_ALL_BRANCHES_1:PMC1",
    /* CLOCK */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "PWR_PKG_ENERGY:PWR0",
    /* DATA */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "MEM_UOPS_RETIRED_LOADS:PMC0,"
    "MEM_UOPS_RETIRED_STORES:PMC1,"
    "UOPS_RETIRED_ALL:PMC2",
    /* ENERGY */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "TEMP_CORE:TMP0,"
    "PWR_PKG_ENERGY:PWR0,"
    "PWR_PP0_ENERGY:PWR1,"
    "PWR_DRAM_ENERGY:PWR3",
    /* ICACHE */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "ICACHE_ACCESSES:PMC0,"
    "ICACHE_MISSES:PMC1,"
    "ICACHE_IFETCH_STALL:PMC2,"
    "ILD_STALL_IQ_FULL:PMC3",
    /* L2 */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "L1D_REPLACEMENT:PMC0,"
    "L2_TRANS_L1D_WB:PMC1,"
    "ICACHE_MISSES:PMC2",
    /* L2CACHE */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "L2_TRANS_ALL_REQUESTS:PMC0,"
    "L2_RQSTS_MISS:PMC1",
    /* L3 */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "L2_LINES_IN_ALL:PMC0,"
    "L2_LINES_OUT_DEMAND_DIRTY:PMC1",
    /* L3CACHE */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "MEM_LOAD_UOPS_RETIRED_L3_ALL:PMC0,"
    "MEM_LOAD_UOPS_RETIRED_L3_MISS:PMC1,"
    "UOPS_RETIRED_ALL:PMC2",
    /* MEM */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "CAS_COUNT_RD:MBOX0C0,"
    "CAS_COUNT_WR:MBOX0C1,"
    "CAS_COUNT_RD:MBOX1C0,"
    "CAS_COUNT_WR:MBOX1C1,"
    "CAS_COUNT_RD:MBOX2C0,"
    "CAS_COUNT_WR:MBOX2C1,"
    "CAS_COUNT_RD:MBOX3C0,"
    "CAS_COUNT_WR:MBOX3C1,"
    "CAS_COUNT_RD:MBOX4C0,"
    "CAS_COUNT_WR:MBOX4C1,"
    "CAS_COUNT_RD:MBOX5C0,"
    "CAS_COUNT_WR:MBOX5C1,"
    "CAS_COUNT_RD:MBOX6C0,"
    "CAS_COUNT_WR:MBOX6C1,"
    "CAS_COUNT_RD:MBOX7C0,"
    "CAS_COUNT_WR:MBOX7C1",
    /* NUMA */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "OFFCORE_RESPONSE_0_LOCAL_DRAM:PMC0,"
    "OFFCORE_RESPONSE_1_REMOTE_DRAM:PMC1",
    /* QPI */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "QBOX_CLOCKTICKS:QBOX0C0,"
    "QBOX_CLOCKTICKS:QBOX1C0,"
    "DIRECT2CORE_SUCCESS_RBT_HIT:QBOX0C1,"
    "DIRECT2CORE_SUCCESS_RBT_HIT:QBOX1C1,"
    "TXL_FLITS_G0_DATA:QBOX0C2,"
    "TXL_FLITS_G0_DATA:QBOX1C2,"
    "TXL_FLITS_G0_NON_DATA:QBOX0C3,"
    "TXL_FLITS_G0_NON_DATA:QBOX1C3",
    /*SBOX */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "RING_BL_USED_ANY:SBOX0C0,"
    "RING_BL_USED_ANY:SBOX1C0,"
    "RING_BL_USED_ANY:SBOX2C0,"
    "RING_BL_USED_ANY:SBOX3C0",
    /* TLB_DATA */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "DTLB_LOAD_MISSES_CAUSES_A_WALK:PMC0,"
    "DTLB_STORE_MISSES_CAUSES_A_WALK:PMC1,"
    "DTLB_LOAD_MISSES_WALK_DURATION:PMC2,"
    "DTLB_STORE_MISSES_WALK_DURATION:PMC3",
    /* TLB_INSTR */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "ITLB_MISSES_CAUSES_A_WALK:PMC0,"
    "ITLB_MISSES_WALK_DURATION:PMC1",
    /* FLOPS_AVX */
    "INSTR_RETIRED_ANY:FIXC0,"
    "CPU_CLK_UNHALTED_CORE:FIXC1,"
    "CPU_CLK_UNHALTED_REF:FIXC2,"
    "AVX_INSTS_CALC:PMC0",
};

const std::vector<const char*> metrics_names[kNumberOfGroups] = {
    /* BRANCH */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Branch rate", "Branch misprediction rate", "Branch misprediction ratio",
     "Instructions per branch"},
    /* CLOCK */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Energy [J]", "Power [W]"},
    /* DATA */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Load to store ratio"},
    /*ENERGY */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Temperature [C]", "Energy [J]", "Power [W]", "Energy PP0 [J]",
     "Power PP0 [W]", "Energy DRAM [J]", "Power DRAM [W]"},
    /* ICACHE */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L1I request rate", "L1I miss rate", "L1I miss ratio", "L1I stalls",
     "L1I stall rate", "L1I queue full stalls", "L1I queue full stall rate"},
    /* L2 */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L2D load bandwidth [MBytes/s]", "L2D load data volume [GBytes]",
     "L2D evict bandwidth [MBytes/s]", "L2D evict data volume [GBytes]",
     "L2 bandwidth [MBytes/s]", "L2 data volume [GBytes]"},
    /* L2CACHE */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L2 request rate", "L2 miss rate", "L2 miss ratio"},
    /* L3 */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L3 load bandwidth [MBytes/s]", "L3 load data volume [GBytes]",
     "L3 evict bandwidth [MBytes/s]", "L3 evict data volume [GBytes]",
     "L3 bandwidth [MBytes/s]", "L3 data volume [GBytes]"},
    /* L3CACHE */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L3 request rate", "L3 miss rate", "L3 miss ratio"},
    /* MEM */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Memory read bandwidth [MBytes/s]", "Memory read data volume [GBytes]",
     "Memory write bandwidth [MBytes/s]", "Memory write data volume [GBytes]",
     "Memory bandwidth [MBytes/s]", "Memory data volume [GBytes]"},
    /* NUMA */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Local DRAM data volume [GByte]", "Local DRAM bandwidth [MByte/s]",
     "Remote DRAM data volume [GByte]", "Remote DRAM bandwidth [MByte/s]",
     "Memory data volume [GByte]", "Memory bandwidth [MByte/s]"},
    /* QPI */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "QPI to LLC data volume [GByte]", "QPI data volume [GByte]",
     "QPI data bandwidth [MByte/s]", "QPI link volume [GByte]",
     "QPI link bandwidth [MByte/s]"},
    /* SBOX */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Ring transfer bandwidth [MByte/s]", "Ring transfer data volume [GByte]"},
    /* TLB_DATA */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L1 DTLB load misses", "L1 DTLB load miss rate",
     "L1 DTLB load miss duration", "L1 DTLB store misses",
     "L1 DTLB store miss rate", "L1 DTLB store miss duration [Cyc]"},
    /* TLB_INSTR */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "L1 ITLB misses", "L1 ITLB miss rate", "L1 ITLB miss duration [Cyc]"},
    /* FLOPS_AVX */
    {"Runtime (RDTSC) [s]", "Runtime unhalted [s]", "Clock [MHz]", "CPI",
     "Packed SP MFLOP/s", "Packed DP MFLOP/s"},
};

const std::vector<std::function<double(int, int)>>
    metrics_functions[kNumberOfGroups] = {
        /* BRANCH */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Branch rate
            [](int group_id, int cpu_id) {
              // PMC0/FIXC0
              return perfmon_getResult(group_id, 3, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Branch misprediction rate
            [](int group_id, int cpu_id) {
              // PMC1/FIXC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Branch misprediction ratio
            [](int group_id, int cpu_id) {
              // PMC1/PMC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 3, cpu_id);
            },
            // Instructions per branch
            [](int group_id, int cpu_id) {
              // FIXC0/PMC0
              return perfmon_getResult(group_id, 0, cpu_id) /
                     perfmon_getResult(group_id, 3, cpu_id);
            },
        },
        /* CLOCK */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Energy [J]
            [](int group_id, int cpu_id) {
              // PWR0
              return perfmon_getResult(group_id, 3, cpu_id);
            },
            // Power [W]
            [](int group_id, int cpu_id) {
              // PWR0/time
              return perfmon_getResult(group_id, 3, cpu_id) /
                     timer_getCpuClock();
            },
        },
        /* DATA */
        {// Runtime (RDTSC) [s]
         [](int group_id, int cpu_id) {
           // time
           return perfmon_getTimeOfGroup(group_id);
         },
         // Runtime unhalted [s]
         [](int group_id, int cpu_id) {
           // FIXC1*inverseClock
           return perfmon_getResult(group_id, 1, cpu_id) / timer_getCpuClock();
           // TODO(guera): likwid 4.1: timer_getCycleClock?
         },
         // Clock [MHz]
         [](int group_id, int cpu_id) {
           // 1.E-06*(FIXC1/FIXC2)/inverseClock
           return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                          perfmon_getResult(group_id, 2, cpu_id)) *
                  timer_getCpuClock();
           // TODO(guera): likwid 4.1: timer_getCycleClock?
         },
         // CPI
         [](int group_id, int cpu_id) {
           // FIXC1/FIXC0
           return perfmon_getResult(group_id, 1, cpu_id) /
                  perfmon_getResult(group_id, 0, cpu_id);

         },
         // Load to store ratio
         [](int group_id, int cpu_id) {
           // PMC0/PMC1
           return perfmon_getResult(group_id, 3, cpu_id) /
                  perfmon_getResult(group_id, 4, cpu_id);
         }},
        /*ENERGY */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Temperature [C]
            [](int group_id, int cpu_id) {
              // TMP0
              return perfmon_getResult(group_id, 3, cpu_id);
            },
            // Energy [J]
            [](int group_id, int cpu_id) {
              // PWR0
              return perfmon_getResult(group_id, 4, cpu_id);
            },
            // Power [W]
            [](int group_id, int cpu_id) {
              // PWR0/time
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getTimeOfGroup(group_id);
            },
            // Energy PP0 [J]
            [](int group_id, int cpu_id) {
              // PWR1
              return perfmon_getResult(group_id, 5, cpu_id);
            },
            // Power PP0 [W]
            [](int group_id, int cpu_id) {
              // PWR1/time
              return perfmon_getResult(group_id, 5, cpu_id) /
                     perfmon_getTimeOfGroup(group_id);
            },
            // Energy DRAM [J]
            [](int group_id, int cpu_id) {
              // PWR3
              return perfmon_getResult(group_id, 6, cpu_id);
            },
            // Power DRAM [W]
            [](int group_id, int cpu_id) {
              // PWR3/time
              return perfmon_getResult(group_id, 6, cpu_id) /
                     perfmon_getTimeOfGroup(group_id);
            },
        },
        /* ICACHE */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1I request rate
            [](int group_id, int cpu_id) {
              // PMC0/FIXC0
              return perfmon_getResult(group_id, 3, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1I miss rate
            [](int group_id, int cpu_id) {
              // PMC1/FIXC0
              return perfmon_getResult(group_id, 5, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1I miss ratio
            [](int group_id, int cpu_id) {
              // PMC1/PMC0
              return perfmon_getResult(group_id, 5, cpu_id) /
                     perfmon_getResult(group_id, 4, cpu_id);
            },
            // L1I stalls
            [](int group_id, int cpu_id) {
              // PMC2
              return perfmon_getResult(group_id, 5, cpu_id);
            },
            // L1I stall rate
            [](int group_id, int cpu_id) {
              // PMC2/FIXC0
              return perfmon_getResult(group_id, 5, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1I queue full stalls
            [](int group_id, int cpu_id) {
              // PMC3
              return perfmon_getResult(group_id, 6, cpu_id);
            },
            // L1I queue full stall rate
            [](int group_id, int cpu_id) {
              // PMC3/FIXC0
              return perfmon_getResult(group_id, 6, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
        },
        /* L2 */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L2D load bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*PMC0*64.0/time
              return 1e-6 * perfmon_getResult(group_id, 3, cpu_id) * 64.0 /
                     perfmon_getTimeOfGroup(group_id);
            },
            // L2D load data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*PMC0*64.0
              return 1e-9 * perfmon_getResult(group_id, 3, cpu_id) * 64.0;
            },
            // L2D evict bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*PMC1*64.0/time
              return 1e-6 * perfmon_getResult(group_id, 4, cpu_id) * 64.0 /
                     perfmon_getTimeOfGroup(group_id);
            },
            // L2D evict data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*PMC1*64.0
              return 1e-9 * perfmon_getResult(group_id, 4, cpu_id) * 64.0;
            },
            // L2 bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*(PMC0+PMC1+PMC2)*64.0/time
              return 1e-6 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id)) *
                     64.0 / perfmon_getTimeOfGroup(group_id);
            },
            // L2 data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*(PMC0+PMC1+PMC2)*64.0
              return 1e-9 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id)) *
                     64.0;
            },
        },
        /* L2CACHE */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L2 request rate
            [](int group_id, int cpu_id) {
              // PMC0/FIXC0
              return perfmon_getResult(group_id, 3, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L2 miss rate
            [](int group_id, int cpu_id) {
              // PMC1/FIXC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L2 miss ratio
            [](int group_id, int cpu_id) {
              // PMC1/PMC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 3, cpu_id);
            },
        },
        /* L3 */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L3 load bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*PMC0*64.0/time
              return 1e-6 * perfmon_getResult(group_id, 3, cpu_id) * 64.0 /
                     perfmon_getTimeOfGroup(group_id);
            },
            // L3 load data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*PMC0*64.0
              return 1e-9 * perfmon_getResult(group_id, 3, cpu_id) * 64.0;
            },
            // L3 evict bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*PMC1*64.0/time
              return 1e-6 * perfmon_getResult(group_id, 4, cpu_id) * 64.0 /
                     perfmon_getTimeOfGroup(group_id);
            },
            // L3 evict data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*PMC1*64.0
              return 1e-9 * perfmon_getResult(group_id, 4, cpu_id) * 64.0;
            },
            // L3 bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*(PMC0+PMC1)*64.0/time
              return 1e-6 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id)) *
                     64.0 / perfmon_getTimeOfGroup(group_id);
            },
            // L3 data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*(PMC0+PMC1)*64.0
              return 1e-9 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id)) *
                     64.0;
            },
        },
        /* L3CACHE */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L3 request rate
            [](int group_id, int cpu_id) {
              // PMC0/PMC2
              return perfmon_getResult(group_id, 3, cpu_id) /
                     perfmon_getResult(group_id, 5, cpu_id);
            },
            // L3 miss rate
            [](int group_id, int cpu_id) {
              // PMC1/PMC2
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 5, cpu_id);
            },
            // L3 miss ratio
            [](int group_id, int cpu_id) {
              // PMC1/PMC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 3, cpu_id);
            },
        },
        /* MEM */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Memory read bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*(MBOX0C0+MBOX1C0+MBOX2C0+
              //          MBOX3C0+MBOX4C0+MBOX5C0+MBOX6C0+MBOX7C0)*64.0/time
              return 1e-6 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 9, cpu_id) +
                             perfmon_getResult(group_id, 11, cpu_id) +
                             perfmon_getResult(group_id, 13, cpu_id) +
                             perfmon_getResult(group_id, 15, cpu_id) +
                             perfmon_getResult(group_id, 17, cpu_id)) *
                     64.0 / perfmon_getTimeOfGroup(group_id);
            },
            // Memory read data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*(MBOX0C0+MBOX1C0+MBOX2C0+MBOX3C0+
              //          MBOX4C0+MBOX5C0+MBOX6C0+MBOX7C0)*64.0
              return 1e-9 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 9, cpu_id) +
                             perfmon_getResult(group_id, 11, cpu_id) +
                             perfmon_getResult(group_id, 13, cpu_id) +
                             perfmon_getResult(group_id, 15, cpu_id) +
                             perfmon_getResult(group_id, 17, cpu_id)) *
                     64.0;
            },
            // Memory write bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*(MBOX0C1+MBOX1C1+MBOX2C1+MBOX3C1+
              //          MBOX4C1+MBOX5C1+MBOX6C1+MBOX7C1)*64.0/time
              return 1e-6 * (perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id) +
                             perfmon_getResult(group_id, 10, cpu_id) +
                             perfmon_getResult(group_id, 12, cpu_id) +
                             perfmon_getResult(group_id, 14, cpu_id) +
                             perfmon_getResult(group_id, 16, cpu_id) +
                             perfmon_getResult(group_id, 18, cpu_id)) *
                     64.0 / perfmon_getTimeOfGroup(group_id);
            },
            // Memory write data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*(MBOX0C1+MBOX1C1+MBOX2C1+MBOX3C1+
              //          MBOX4C1+MBOX5C1+MBOX6C1+MBOX7C1)*64.0
              return 1e-9 * (perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id) +
                             perfmon_getResult(group_id, 10, cpu_id) +
                             perfmon_getResult(group_id, 12, cpu_id) +
                             perfmon_getResult(group_id, 14, cpu_id) +
                             perfmon_getResult(group_id, 16, cpu_id) +
                             perfmon_getResult(group_id, 18, cpu_id)) *
                     64.0;
            },
            // Memory bandwidth [MBytes/s]
            [](int group_id, int cpu_id) {
              // 1.0E-06*(MBOX0C0+MBOX1C0+MBOX2C0+MBOX3C0+
              //          MBOX4C0+MBOX5C0+MBOX6C0+MBOX7C0+
              //          MBOX0C1+MBOX1C1+MBOX2C1+MBOX3C1+
              //          MBOX4C1+MBOX5C1+MBOX6C1+MBOX7C1)*64.0/time
              return 1e-6 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id) +
                             perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id) +
                             perfmon_getResult(group_id, 9, cpu_id) +
                             perfmon_getResult(group_id, 10, cpu_id) +
                             perfmon_getResult(group_id, 11, cpu_id) +
                             perfmon_getResult(group_id, 12, cpu_id) +
                             perfmon_getResult(group_id, 13, cpu_id) +
                             perfmon_getResult(group_id, 14, cpu_id) +
                             perfmon_getResult(group_id, 15, cpu_id) +
                             perfmon_getResult(group_id, 16, cpu_id) +
                             perfmon_getResult(group_id, 17, cpu_id) +
                             perfmon_getResult(group_id, 18, cpu_id)) *
                     64.0 / perfmon_getTimeOfGroup(group_id);
            },
            // Memory data volume [GBytes]
            [](int group_id, int cpu_id) {
              // 1.0E-09*(MBOX0C0+MBOX1C0+MBOX2C0+MBOX3C0+
              //          MBOX4C0+MBOX5C0+MBOX6C0+MBOX7C0+
              //          MBOX0C1+MBOX1C1+MBOX2C1+MBOX3C1+
              //          MBOX4C1+MBOX5C1+MBOX6C1+MBOX7C1)*64.0
              return 1e-9 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id) +
                             perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id) +
                             perfmon_getResult(group_id, 9, cpu_id) +
                             perfmon_getResult(group_id, 10, cpu_id) +
                             perfmon_getResult(group_id, 11, cpu_id) +
                             perfmon_getResult(group_id, 12, cpu_id) +
                             perfmon_getResult(group_id, 13, cpu_id) +
                             perfmon_getResult(group_id, 14, cpu_id) +
                             perfmon_getResult(group_id, 15, cpu_id) +
                             perfmon_getResult(group_id, 16, cpu_id) +
                             perfmon_getResult(group_id, 17, cpu_id) +
                             perfmon_getResult(group_id, 18, cpu_id)) *
                     64.0;
            },
        },
        /* NUMA */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Local DRAM data volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-09*PMC0*64
              return 1e-9 * perfmon_getResult(group_id, 3, cpu_id) * 64.0;
            },
            // Local DRAM bandwidth [MByte/s]
            [](int group_id, int cpu_id) {
              // 1.E-06*(PMC0*64)/time
              return 1e-6 * (perfmon_getResult(group_id, 3, cpu_id)) /
                     perfmon_getTimeOfGroup(group_id);
            },
            // Remote DRAM data volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-09*PMC1*64
              return 1e-9 * perfmon_getResult(group_id, 4, cpu_id) * 64.0;
            },
            // Remote DRAM bandwidth [MByte/s]
            [](int group_id, int cpu_id) {
              // 1.E-06*(PMC1*64)/time
              return 1e-6 * (perfmon_getResult(group_id, 4, cpu_id) * 64.0) /
                     perfmon_getTimeOfGroup(group_id);
            },
            // Memory data volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-09*(PMC0+PMC1)*64
              return 1e-9 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id)) *
                     64.0;
            },
            // Memory bandwidth [MByte/s]
            [](int group_id, int cpu_id) {
              // 1.E-06*((PMC0+PMC1)*64)/time
              return 1e-6 * ((perfmon_getResult(group_id, 3, cpu_id) +
                              perfmon_getResult(group_id, 4, cpu_id)) *
                             64.0) /
                     perfmon_getTimeOfGroup(group_id);
            },
        },
        /* QPI */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // QPI to LLC data volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-09*(QBOX0C1+QBOX1C1)*64
              return 1e-9 * (perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id)) *
                     64.0;
            },
            // QPI data volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-06*(QBOX0C2+QBOX1C2)*8
              return 1e-6 * (perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id)) *
                     8.0;
            },
            // QPI data bandwidth [MByte/s]
            [](int group_id, int cpu_id) {
              // 1.E-09*(QBOX0C2+QBOX1C2)*8/time
              return 1e-9 * (perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id)) *
                     8.0 / perfmon_getTimeOfGroup(group_id);
            },
            // QPI link volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-06*(QBOX0C2+QBOX1C2+QBOX0C3+QBOX1C3)*8
              return 1e-6 * (perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id) +
                             perfmon_getResult(group_id, 9, cpu_id) +
                             perfmon_getResult(group_id, 10, cpu_id)) *
                     8.0;
            },
            // QPI link bandwidth [MByte/s]
            [](int group_id, int cpu_id) {
              // 1.E-09*(QBOX0C2+QBOX1C2+QBOX0C3+QBOX1C3)*8/time
              return 1e-9 * (perfmon_getResult(group_id, 7, cpu_id) +
                             perfmon_getResult(group_id, 8, cpu_id) +
                             perfmon_getResult(group_id, 9, cpu_id) +
                             perfmon_getResult(group_id, 10, cpu_id)) *
                     8.0 / perfmon_getTimeOfGroup(group_id);
            },
        },
        /* SBOX */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Ring transfer bandwidth [MByte/s]
            [](int group_id, int cpu_id) {
              // 1.E-06*(SBOX0C0+SBOX1C0+SBOX2C0+SBOX3C0)*32/time
              return 1e-6 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id)) *
                     32.0 / perfmon_getTimeOfGroup(group_id);
            },
            // Ring transfer data volume [GByte]
            [](int group_id, int cpu_id) {
              // 1.E-09*(SBOX0C0+SBOX1C0+SBOX2C0+SBOX3C0)*32
              return 1e-9 * (perfmon_getResult(group_id, 3, cpu_id) +
                             perfmon_getResult(group_id, 4, cpu_id) +
                             perfmon_getResult(group_id, 5, cpu_id) +
                             perfmon_getResult(group_id, 6, cpu_id)) *
                     32.0;
            },
        },
        /* TLB_DATA */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1 DTLB load misses
            [](int group_id, int cpu_id) {
              // PMC0
              return perfmon_getResult(group_id, 3, cpu_id);
            },
            // L1 DTLB load miss rate
            [](int group_id, int cpu_id) {
              // PMC0/FIXC0
              return perfmon_getResult(group_id, 3, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1 DTLB load miss duration
            [](int group_id, int cpu_id) {
              // PMC2/PMC0
              return perfmon_getResult(group_id, 5, cpu_id) /
                     perfmon_getResult(group_id, 3, cpu_id);
            },
            // L1 DTLB store misses
            [](int group_id, int cpu_id) {
              // PMC1
              return perfmon_getResult(group_id, 4, cpu_id);
            },
            // L1 DTLB store miss rate
            [](int group_id, int cpu_id) {
              // PMC1/FIXC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1 DTLB store miss duration [Cyc]
            [](int group_id, int cpu_id) {
              // PMC3/PMC1
              return perfmon_getResult(group_id, 6, cpu_id) /
                     perfmon_getResult(group_id, 4, cpu_id);
            },
        },
        /* TLB_INSTR */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1 ITLB misses
            [](int group_id, int cpu_id) {
              // PMC0
              return perfmon_getResult(group_id, 3, cpu_id);
            },
            // L1 ITLB miss rate
            [](int group_id, int cpu_id) {
              // PMC0/FIXC0
              return perfmon_getResult(group_id, 3, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // L1 ITLB miss duration [Cyc]
            [](int group_id, int cpu_id) {
              // PMC1/PMC0
              return perfmon_getResult(group_id, 4, cpu_id) /
                     perfmon_getResult(group_id, 3, cpu_id);
            },
        },
        /* FLOPS_AVX */
        {
            // Runtime (RDTSC) [s]
            [](int group_id, int cpu_id) {
              // time
              return perfmon_getTimeOfGroup(group_id);
            },
            // Runtime unhalted [s]
            [](int group_id, int cpu_id) {
              // FIXC1*inverseClock
              return perfmon_getResult(group_id, 1, cpu_id) /
                     timer_getCpuClock();
            },
            // Clock [MHz]
            [](int group_id, int cpu_id) {
              // 1.E-06*(FIXC1/FIXC2)/inverseClock
              return 1e-6 * (perfmon_getResult(group_id, 1, cpu_id) /
                             perfmon_getResult(group_id, 2, cpu_id)) *
                     timer_getCpuClock();
            },
            // CPI
            [](int group_id, int cpu_id) {
              // FIXC1/FIXC0
              return perfmon_getResult(group_id, 1, cpu_id) /
                     perfmon_getResult(group_id, 0, cpu_id);
            },
            // Packed SP MFLOP/s
            [](int group_id, int cpu_id) {
              // 1.0E-06*(PMC0*8.0)/time
              return 1e-6 * perfmon_getResult(group_id, 3, cpu_id) * 8 /
                     perfmon_getTimeOfGroup(group_id);
            },
            // Packed DP MFLOP/s
            [](int group_id, int cpu_id) {
              // 1.0E-06*(PMC0*8.0)/time
              return 1e-6 * perfmon_getResult(group_id, 3, cpu_id) * 4 /
                     perfmon_getTimeOfGroup(group_id);
            },
        },
};

void test() {
  const int event = 1;
  int handle = perfmon_addEventSet(const_cast<char*>(eventsets[event]));
  perfmon_setupCounters(handle);
  perfmon_startCounters();
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  perfmon_startCounters();
  std::cout << metrics_functions[event][4](handle, 0) << std::endl;

  // perfmon_setupCounters(handle);
  perfmon_startCounters();
  std::this_thread::sleep_for(std::chrono::milliseconds(2000));
  perfmon_startCounters();
  std::cout << metrics_functions[event][4](handle, 0) << std::endl;

  perfmon_setupCounters(handle);
  perfmon_startCounters();
  std::this_thread::sleep_for(std::chrono::milliseconds(3000));
  perfmon_startCounters();
  std::cout << metrics_functions[event][4](handle, 0) << std::endl;
}

}  // namespace

namespace exahype {
namespace profilers {
namespace likwid {

LikwidPerformanceMonitoringModule::LikwidPerformanceMonitoringModule(
    const LikwidProfilerState& state, const std::string& group_name)
    : LikwidModule(state), group_name_(group_name) {
  int errcode;

  // Initialize perfmon module
  errcode = perfmon_init(1, const_cast<int*>(&state_.cpu_));
  if (errcode != 0) {
    std::cerr << "LikwidPerformanceMonitoringModule: perfmon_init returned "
              << errcode << " != 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  group_index_ = std::distance(
      std::begin(groups), std::find_if(std::begin(groups), std::end(groups),
                                       [this](const char* group) {
                                         return group_name_.compare(group) == 0;
                                       }));
  if (group_index_ == std::distance(std::begin(groups), std::end(groups))) {
    std::cerr << "LikwidPerformanceMonitoringModule: group_name_ = "
              << group_name_ << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  timer_init();

  // test();
  // std::exit(0);
}

LikwidPerformanceMonitoringModule::~LikwidPerformanceMonitoringModule() {
  perfmon_finalize();
  timer_finalize();
}

void LikwidPerformanceMonitoringModule::setNumberOfTags(int n) {
  group_handles_.reserve(n);
  counts_.reserve(n);
}

void LikwidPerformanceMonitoringModule::registerTag(const std::string& tag) {
  int handle = perfmon_addEventSet(const_cast<char*>(eventsets[group_index_]));
  if (handle < 0) {
    std::cerr << "LikwidPerformanceMonitoringModule: addEventSet returned "
              << handle << " < 0 for '" << eventsets[group_index_] << "'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  group_handles_[tag] = handle;
  counts_[tag] = 0;
}

void LikwidPerformanceMonitoringModule::start(const std::string& tag) {
  int errcode;

  errcode = perfmon_setupCounters(group_handles_.at(tag));
  assert((errcode >= 0) &&
         "LikwidPerformanceMonitoringModule: setupCounters failed");

  errcode = perfmon_startCounters();
  assert((errcode == 0) &&
         "LikwidPerformanceMonitoringModule: startCounters failed");
}

void LikwidPerformanceMonitoringModule::stop(const std::string& tag) {
  // Note: Tagged regions may not overlap, since stopCounters stops all groups
  int errcode;

  errcode = perfmon_stopCounters();
  assert((errcode == 0) &&
         "LikwidPerformanceMonitoringModule: stopCounters failed");

  counts_[tag]++;
}

void LikwidPerformanceMonitoringModule::writeToOstream(std::ostream* os) const {
  assert((metrics_names[group_index_].size() ==
          metrics_functions[group_index_].size()) &&
         "LikwidPerformanceMonitoringModule: metrics_names and "
         "metrics_functions don't have the same length");

  // For all tags
  for (const auto& tag_group_handle_pair : group_handles_) {
    *os << "PerformanceMonitoringModule: " << tag_group_handle_pair.first
        << " count " << counts_.at(tag_group_handle_pair.first) << std::endl;

    int number_of_metrics = metrics_names[group_index_].size();
    // for all metrics
    for (int i = 0; i < number_of_metrics; i++) {
      *os << "PerformanceMonitoringModule: " << tag_group_handle_pair.first
          << " " << metrics_names[group_index_][i] << " "
          << metrics_functions[group_index_]
                              [i](tag_group_handle_pair.second, state_.cpu_)
          << std::endl;
    }

    *os << "PerformanceMonitoringModule: " << tag_group_handle_pair.first
        << " runtime / count "
        << metrics_functions[group_index_]
                            [0](tag_group_handle_pair.second, state_.cpu_) /
               counts_.at(tag_group_handle_pair.first)
        << std::endl;

    *os << "PerformanceMonitoringModule: " << tag_group_handle_pair.first
        << " runtime unhalted / count "
        << metrics_functions[group_index_]
                            [1](tag_group_handle_pair.second, state_.cpu_) /
               counts_.at(tag_group_handle_pair.first)
        << std::endl;

    std::istringstream iss(eventsets[group_index_]);
    std::string counter;
    int counter_index = 0;

    // For all counters
    while (getline(iss, counter, ',')) {
      *os << "PerformanceMonitoringModule: " << tag_group_handle_pair.first
          << " " << counter << " "
          << perfmon_getResult(tag_group_handle_pair.second, counter_index,
                               state_.cpu_)
          << std::endl;
      counter_index++;
    }
  }
}

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE
