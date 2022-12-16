// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <frost/frost_parameter_interface.hpp>
#include <frost/sampler.hpp>

namespace frost {

void build_trimesh3( frantic::geometry::trimesh3& outMesh, std::size_t& outParticleCount,
                     frost_parameter_interface& params,
                     boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                     frantic::logging::progress_logger& progressLogger );

void populate_mesh_channels( frantic::geometry::trimesh3& outMesh, frost_parameter_interface& params,
                             boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                             const frantic::channels::channel_propagation_policy& cpp,
                             frantic::logging::progress_logger& progressLogger );

void evolve_mesh( frantic::geometry::trimesh3& outMesh, frost_parameter_interface& params, int iterations,
                  float relativeSpacing, float relaxWeight,
                  boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                  const frantic::channels::channel_propagation_policy& cpp,
                  frantic::logging::progress_logger& progressLogger );

frost::sampler::ptr_type create_sampler( frost_parameter_interface& params,
                                         boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                         frantic::logging::progress_logger& progressLogger );

} // namespace frost
