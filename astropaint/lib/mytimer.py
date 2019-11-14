#!/usr/bin/env python

__author__ = "Siavash Yasini"

from contextlib import contextmanager

@contextmanager
def timeit(process_name="Process"):
    """Time the code in mins"""
    import time

    time_stamp = time.strftime("%H:%M:%S %p")
    print("{:=>50}\n{} started at {}\n".format("", process_name, time_stamp))
    t_i = time.time()

    yield

    t_f = time.time()

    t = t_f-t_i

    print("\n{} was done in {:.1f} min.\n{:=>50}\n".format(process_name, t/60, ""))