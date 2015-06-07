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
"""
import os
import errno


def makedirs(name, mode=0777):
    """
    Create the directory if it doesn't already exist

    NB: mimics the os.makedirs function

    :param directory_path:
    :param directory_path:
    :return:
    """
    try:
        os.makedirs(name, mode)
    except OSError as err:
        # be happy if someone already created the path
        if err.errno == errno.EEXIST and os.path.isdir(name):
            pass
        else:
            raise

