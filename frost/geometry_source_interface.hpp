// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/trimesh3.hpp>

class geometry_source_interface {
  public:
    virtual bool get_geometry( frantic::geometry::trimesh3& outMesh, double timeInSeconds,
                               std::size_t shapeNumber ) = 0;
    virtual std::size_t get_geometry_count() = 0;
};
