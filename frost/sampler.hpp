// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/shared_ptr.hpp>

#include <frantic/channels/channel_map.hpp>

namespace frost {

class sampler {
  public:
    typedef boost::shared_ptr<sampler> ptr_type;

    virtual ~sampler() {}

    virtual const frantic::channels::channel_map& get_channel_map() const = 0;

    virtual void get_particles_in_range( const frantic::graphics::vector3f& position,
                                         std::vector<char*>& outParticles ) = 0;

    virtual void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                                     std::vector<float>& outWeights ) = 0;
};

} // namespace frost
