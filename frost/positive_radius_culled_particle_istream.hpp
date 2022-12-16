// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/none.hpp>

#include <frantic/particles/streams/culling_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

class positive_radius_culling_policy
    : public frantic::particles::streams::culling_policy_base<positive_radius_culling_policy> {
    frantic::channels::channel_cvt_accessor<float> m_radiusAccessor;

  public:
    typedef boost::none_t args_type;

    /**
     * Constructs an id culling policy, whose arguments are a set of ids, and a boolean of whether to cull
     * particles with ids that belong to the set, or do not belong to the set
     *
     * @param args tuple of the policy arguments for the policy
     * @param pcm the channel map to use for the stream
     */
    positive_radius_culling_policy( const args_type& /*args*/, const frantic::channels::channel_map& pcm ) {
        set_channel_map( pcm );
    }

    // Need a splitting constructor to support TBB
    positive_radius_culling_policy( const positive_radius_culling_policy& lhs, tbb::split )
        : m_radiusAccessor( lhs.m_radiusAccessor ) {}

    /**
     * Sets the channel map with which the particles will be read from the stream
     *
     * @param pcm a channel map
     */
    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_radiusAccessor = pcm.get_cvt_accessor<float>( _T("Radius") );
    }

    /**
     * @param particle - a data pointer to the particle
     * @returns true if the particle should be culled
     */
    bool cull( char* particle ) const { return m_radiusAccessor.get( particle ) <= 0; }
};

typedef frantic::particles::streams::culling_particle_istream<positive_radius_culling_policy>
    positive_radius_culled_particle_istream;
