#pragma once
#include "common.h"

struct PolymerBase
{
    int output_flag = 0;
    int length = 0;
};

struct SequenceSegment
{
    SpeciesID monomer;
    size_t length;
};

struct SegmentData
{
    size_t count = 0;
    size_t segments = 0;
    size_t segment_lengths_squared = 0;

    void add_segment(SequenceSegment segment)
    {
        count += segment.length;
        ++segments;
        segment_lengths_squared += segment.length * segment.length;
    }

    void remove_segment(SequenceSegment segment)
    {
        count -= segment.length;
        --segments;
        segment_lengths_squared -= segment.length * segment.length;
    }
};

template <int _NUM_MONOMERS>
struct Fragment
{
    
};

template <int _NUM_MONOMERS>
struct Microstructure
{
    std::vector<Fragment> fragments;
};

template <int _NUM_MONOMERS>
struct CopolymerBase : public PolymerBase
{
    int counts[_NUM_MONOMERS];
    int segments[_NUM_MONOMERS];
    int segment_lengths_squared[_NUM_MONOMERS];
};
