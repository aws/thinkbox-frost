// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

class create_id_channel_from_index_particle_istream : public frantic::particles::streams::delegated_particle_istream {
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_idAccessor;

    frantic::tstring m_idChannelName;

    frantic::channels::channel_map m_nativeMap; // May have to add to the native map;

  public:
    create_id_channel_from_index_particle_istream( boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                                   const frantic::tstring& idChannelName )
        : delegated_particle_istream( pin )
        , m_idChannelName( idChannelName ) {
        set_channel_map( m_delegate->get_channel_map() );

        m_nativeMap = m_delegate->get_native_channel_map();
        if( !m_nativeMap.has_channel( idChannelName ) ) {
            m_nativeMap.append_channel<boost::int32_t>( idChannelName );
        }
    }

    virtual ~create_id_channel_from_index_particle_istream() {}

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_delegate->set_channel_map( pcm );
        m_idAccessor.reset();
        if( pcm.has_channel( m_idChannelName ) ) {
            m_idAccessor = pcm.get_cvt_accessor<boost::int32_t>( m_idChannelName );
        }
    }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    bool get_particle( char* p ) {
        if( !m_delegate->get_particle( p ) ) {
            return false;
        }

        const boost::int64_t particleIndex = m_delegate->particle_index();

        m_idAccessor.set( p, static_cast<boost::int32_t>( particleIndex ) );

        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        const boost::int64_t beforeIndex = m_delegate->particle_index();

        const bool notEos = m_delegate->get_particles( buffer, numParticles );

        const std::size_t structureSize = m_delegate->get_channel_map().structure_size();

        boost::int64_t particleIndex = beforeIndex;
        for( std::size_t i = 0; i < numParticles; ++i ) {
            ++particleIndex;
            m_idAccessor.set( buffer, static_cast<boost::int32_t>( particleIndex ) );
            buffer += structureSize;
        }

        return notEos;
    }
};
