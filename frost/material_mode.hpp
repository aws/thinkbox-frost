// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frost {

namespace material_mode {
enum option {
    single,
    mtlindex_channel,
    shape_number,
    material_id_from_geometry,
    material_from_geometry,
    //
    count
};
}

} // namespace frost
