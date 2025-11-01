// /*

// Brainstorming:

// `Polymer` class is a representation of a macromolecule made up of monomers.

// It can be linear, branched, crosslinked, etc.

// `Segment` struct is a linear, contiguous block of identical monomers within a copolymer.
// `SegmentStats` struct accumulates statistics about segments in a copolymer.

// `Fragment` struct represents a linear section of polymer

// `Polymer` is templated from:
// - Different polymer architecture (linear, branched, ...)
// - Different segment storage strategies (full sequence, compressed, ...)
// - - Full sequence (vector of SpeciesID)
// - - Sequence Statistics (SegmentStats)



// */

// #pragma once

// #include "common.h"

// struct PolymerBase
// {
//     int output_flag = 0;
//     int length = 0;
// };

// struct Homopolymer : public PolymerBase
// {
//     PolymerBase to_data() const { return {output_flag, length}; }
// };

// template <int _NUM_MONOMERS>
// struct CopolymerBase : public PolymerBase
// {
//     std::vector<uint32_t> monomer_counts;
// };

// template <int _NUM_MONOMERS>
// struct SequenceCopolymer : CopolymerBase<_NUM_MONOMERS>
// {
//     std::vector<uint32_t> segments;
//     std::vector<uint32_t> segment_lengths_squared;
// };

// template <int _NUM_MONOMERS>
// struct Polymer
// {
//     using ChainType = void;
// };

// template <>
// struct Polymer<1>
// {
//     using ChainType = Homopolymer;
// };

// template <>
// struct Polymer<2>
// {
//     using ChainType = SequenceCopolymer<2>;
// };

// struct Segment
// {
//     SpeciesID monomer;
//     uint32_t length;
// };


// struct SegmentStats
// {
//     size_t count = 0;
//     size_t segments = 0;
//     size_t segment_lengths_squared = 0;

//     void add_segment(Segment segment)
//     {
//         count += segment.length;
//         ++segments;
//         segment_lengths_squared += segment.length * segment.length;
//     }

//     void remove_segment(Segment segment)
//     {
//         count -= segment.length;
//         --segments;
//         segment_lengths_squared -= segment.length * segment.length;
//     }
// };

// template <int _NUM_MONOMERS>
// struct CopolymerData
// {
//     int output_flag = 0;
//     int length = 0;
//     int counts[_NUM_MONOMERS];
//     int segments[_NUM_MONOMERS];
//     int segment_lengths_squared[_NUM_MONOMERS];
// };

// enum class ConnectType
// {
//     SCB,
//     LCB,
// };

// template <int _NUM_MONOMERS>
// struct Connect
// {
//     Polymer<_NUM_MONOMERS> *chain1;
//     Polymer<_NUM_MONOMERS> *chain2;
//     ConnectType type;
// };

// template <int _NUM_MONOMERS>
// class Fragment
// {
//     using chain = void;
// };

// template <>
// class Fragment<1>
// {
//     Homopolymer chain;

//     void add_monomer(SpeciesID monomer) { chain.length++; }
//     void remove_monomer()
//     {
//         if (chain.length > 0)
//             chain.length--;
//     }
// };

// template <>
// class Fragment<2>
// {
//     CopolymerBase chain;

//     void add_monomer(SpeciesID monomer) { chain.add_monomer(monomer); }
//     void remove_monomer() { chain.remove_monomer(); }
// };

// template <int SIZE>
// struct SegmentBuffer
// {
//     Segment segments[SIZE];
//     int count = 0;

//     void push(const Segment &segment)
//     {
//         if (count >= SIZE)
//             console::error("SegmentBuffer overflow.");
//         segments[count++] = segment;
//     }

//     Segment pop()
//     {
//         if (count <= 0)
//             return {INVALID_SPECIES_ID, 0};
//         return segments[--count];
//     }

//     bool empty() const { return count == 0; }
// };

// class CopolymerBase : public PolymerBase
// {
//     Segment current_segment;
//     SegmentBuffer<10> buffer;

// public:
//     void add_monomer(SpeciesID monomer)
//     {
//         if (monomer == INVALID_SPECIES_ID)
//             return;

//         this->length++;
//         if (current_segment.monomer == monomer)
//             current_segment.length++;
//         else
//         {
//             if (current_segment.length > 0)
//                 buffer.push(current_segment);
//             current_segment = {monomer, 1};
//         }
//     }

//     SpeciesID remove_monomer()
//     {
//         if (length == 0)
//             return INVALID_SPECIES_ID;

//         this->length--;

//         if (current_segment.length > 1)
//             current_segment.length--;
//         else
//             current_segment = buffer.pop();

//         return current_segment.monomer;
//     }
// };

// template <int NUM_MONOMERS>
// struct PolymerSpecies
// {
//     using Polymer = void;
// };

// template <>
// struct PolymerSpecies<1>
// {
//     using ChainType = Homopolymer;
//     struct Polymer
//     {
//         ChainType chain;
//     };
// };

// struct ChainView
// {
//     uint32_t length() const;
//     std::vector<const Segment> segments() const;
// };