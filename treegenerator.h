#ifndef STUDENT_H
#define STUDENT_H

#include <vector>

enum node_type {AXILLARY_BUD, TERMINAL_BUD, BRANCH, BASE};

typedef float vec3_t[3];

typedef struct node
{
	node_type type;
	vec3_t location;
	vec3_t dir;
	bool has_branch;
	node *branch_node;
	double light_info;
	double avg_light_info;
	int total_buds;
	double base_v;
} node;

typedef struct branch
{
	std::vector<node *> nodes;
	std::vector<node *> sorted_nodes;
	node *base_node;
} base_node;


class TreeGenerator
{

private:
	float ***shadow_grid;
	int q_max;
	float a, c;
	int b;
	float internode_length;
	
	void min_voxel(int *, std::vector<int *>);
	double vector_angle(vec3_t, vec3_t);
	void voxel_from_location(int *, vec3_t);
	double weight_function(int, double, double, double, int);
	double get_resources(branch, node*);
	
	void grow_metamors(branch *, node *, int, double);
	void locations_from_voxel(int *, vec3_t *);
	void simple_location_from_voxel(int *, vec3_t);
	void calculate_location(node *, node *, double);
	void normalize(vec3_t);
	void init_tree();
	void grow_tree_step();
	void calculate_growth_direction(node *, double);
	void get_bud_direction(vec3_t, vec3_t);
public:
	//4-1
	TreeGenerator();
	std::vector<branch> branches;
	void init_shadow_grid();
	void propogate_shadow(int, int, int);
	float get_light(int, int, int);
	void calculate_optimal_direction(node *, double, vec3_t);
	void free_grid();
	//4-2
	void priority_model_step();
	void propogate_light_info();
	void grow_tree(int);
	void print_tree();

};

#endif