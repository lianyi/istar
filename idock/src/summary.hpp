#pragma once
#ifndef IDOCK_SUMMARY_HPP
#define IDOCK_SUMMARY_HPP

#include "conformation.hpp"

/// Represents a summary of docking results of a ligand.
class summary
{
public:
	size_t index;
	fl energy;
	fl rfscore;
	conformation conf;
	explicit summary(const size_t index, const fl energy, const fl rfscore, const conformation& conf) : index(index), energy(energy), rfscore(rfscore), conf(conf)
	{
	}

	summary(const summary&) = default;
	summary(summary&&) = default;
	summary& operator=(const summary&) = default;
	summary& operator=(summary&&) = default;
};

/// For sorting ptr_vector<summary>.
inline bool operator<(const summary& a, const summary& b)
{
	return a.energy < b.energy;
//	return a.rfscore > b.rfscore;
}

#endif
