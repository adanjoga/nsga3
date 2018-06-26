#ifndef REFPOINT_H_
#define REFPOINT_H_

#include "global.h"

class TRefpoint {
private:
	size_t num_members;			// number of members associated to the point
	vector<pair<size_t, double> > candidate_members; // candidate members

public:
	TRefpoint();
	TRefpoint(size_t M);
	virtual ~TRefpoint();

	vector<double> refpoint; 	// consider to move this to private section

	size_t member_size();
	bool has_candidate_members();
	int get_nearest_member();
	int get_random_member();

	void sum_member();
	void add_member(size_t idx, double distance);
	void remove_member(size_t idx);
	void clear_refpoint();

	void show_rpoint();
};

TRefpoint::TRefpoint() {
	num_members = 0;
}

TRefpoint::TRefpoint(size_t M) {
	num_members = 0;
	for (size_t i = 0; i < M; ++i)
		refpoint.push_back(0.0);
}

TRefpoint::~TRefpoint() {}

size_t TRefpoint::member_size() {
	return num_members;
}

bool TRefpoint::has_candidate_members() {
	return !candidate_members.empty();
}

int TRefpoint::get_nearest_member() {

	int min_indv = -1;
	double min_dist = MAX_DOUBLE;
	for (size_t i = 0; i < candidate_members.size(); ++i) {
		if (candidate_members[i].second < min_dist) {
			min_dist = candidate_members[i].second;
			min_indv = candidate_members[i].first;
		}
	}
	return min_indv;
}

int TRefpoint::get_random_member() {
	return candidate_members[rand() % candidate_members.size()].first;
}

void TRefpoint::sum_member() {
	num_members += 1;
}

void TRefpoint::add_member(size_t idx, double distance) {
	candidate_members.push_back(make_pair(idx, distance));
}

void TRefpoint::remove_member(size_t idx) {
	for (size_t i = 0; i < candidate_members.size(); ++i) {
		if (candidate_members[i].first == idx) {
			candidate_members.erase(candidate_members.begin() + i);
			return;
		}
	}
}

void TRefpoint::clear_refpoint() {
	num_members = 0;
	candidate_members.clear();
}

void TRefpoint::show_rpoint() {
	for (size_t i = 0; i < numObjectives; i++)
		printf("%f ", refpoint[i]);
	printf("\n");
}


#endif /* REFPOINT_H_ */
