// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frost/meshing_policy.hpp>

class frost_parameter_interface;

meshing_policy::ptr_type create_meshing_policy( const frost_parameter_interface& params );
