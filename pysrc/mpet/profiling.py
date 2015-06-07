"""!
@file
@date 07 Jun 2015

@license
Copyright 2015 Brett Tully

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@description
A wrapper around the google-perftools libraries
"""
import ctypes.util
import os
import mpet.io

try:
    cpu_profiler = ctypes.CDLL(ctypes.util.find_library('profiler'))
except ImportError as err:
    print err
    raise ImportError("Could not load the CPU profiler, have you install googler-perftools?")

try:
    memory_profiler = ctypes.CDLL(ctypes.util.find_library('tcmalloc'))
except ImportError as err:
    print err
    raise ImportError("Could not load the memory profiler, have you install googler-perftools?")


def start(output_prefix="mybin", profile_memory=False, cpu_profile_freq=1000):
    """
    Start the profiling tools

    @param output_prefix optional string to describe the prefix to the profiling files
    """
    os.environ["CPUPROFILE_FREQUENCY"] = str(int(cpu_profile_freq))
    output_path = os.path.dirname(output_prefix)
    if output_path is not '':
        mpet.io.makedirs(output_path)
    result = cpu_profiler.ProfilerStart(output_prefix + ".prof")
    if result < 0:
        raise RuntimeError("CPU profiling failed to start: %s" % result)

    if profile_memory:
        memory_profiler.HeapProfilerStart(output_prefix)
        if not memory_profiler.IsHeapProfilerRunning():
            raise RuntimeError("Memory profiling failed to start")


def profiling_memory_on():
    """
    Returns true if memory profiling is switched on
    """
    return memory_profiler.IsHeapProfilerRunning()


def stop():
    """
    Stop the profiling tools
    """
    cpu_profiler.ProfilerStop()
    if profiling_memory_on():
        memory_profiler.HeapProfilerDump()
        memory_profiler.HeapProfilerStop()

