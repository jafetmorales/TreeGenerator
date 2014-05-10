#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "treegenerator.h"





TreeGenerator::TreeGenerator()
{
	a = 0.2;
	c = 1.0;
	b = 2;
	q_max = 4;
	srand(time(NULL));
}	

void TreeGenerator::init_shadow_grid()
{

	//initialize a 200 x 200 x 200 shadow grid with each voxel initialized to 0.0
	shadow_grid = new float**[200];
	int i, j, k;
	for(i=0; i<200; i++)
	{
		shadow_grid[i] = new float*[200];
		for(j=0; j<200; j++)
		{
			shadow_grid[i][j] = new float[200];
			for(k=0; k<200; k++)
			{
				shadow_grid[i][j][k] = 0.0;
			}
		}
	}
	//printf("Shadow grid initialized\n");

}

void TreeGenerator::init_tree()
{
	//printf("Initializing Tree\n");
	//make base and terminal node - both same
	node *base = new node;
	node *start = new node;
	base->type = BASE;
	start->type = TERMINAL_BUD;
	base->location[0] = start->location[0] = base->location[2] = start->location[2] = 0.0;
	base->location[1] = 0.0;
	start->location[1] = 1.0;
	base->dir[0] = start->dir[0] = base->dir[2] = start->dir[2] = 0.0;
	base->dir[1] = start->dir[1] = 1.0; //stright up
	start->has_branch = false;

	//add them to a new branch
	branch *first_branch = new branch;
	first_branch->nodes.push_back(base);

	first_branch->nodes.push_back(start);

	branches.push_back(*first_branch);

	//make sure we get the shadow (not that it matters a whole lot)
	propogate_shadow(99,1,99);

	int bv[3], sv[3];
	voxel_from_location(bv, branches.at(0).nodes.at(0)->location);
	voxel_from_location(sv, branches.at(0).nodes.at(1)->location);

	//printf("Base location: (%d, %d, %d)\n", bv[0], bv[1], bv[2]);
	//printf("Start location: (%d, %d, %d)\n", sv[0], sv[1], sv[2]);

}

void TreeGenerator::propogate_shadow(int I, int J, int K)
{
	//printf("propogating shadow\n");
	//printf("Start voxel = [%d][%d][%d]\n", I, J, K);
	for(int q=0; q<q_max; q++)
	{	
		if((J-q) < 0) break; //out of bounds
		float delta_s = a * pow(b, (double) -q);
		for(int pi = -q; pi <= q; pi++)
		{
			for(int pk = -q; pk <= q; pk++)
			{
				int curr_i = I+pi;
				int curr_j = J-q;
				int curr_k = K+pk;
				if( curr_i < 0 || curr_i >= 200 ||
					curr_j < 0 || curr_j >= 200 ||
					curr_k < 0 || curr_k >= 200 )
						continue;
				//printf("propogate_shadow:\tHitting voxel ( %d, %d, %d )\n", curr_i, curr_j, curr_k);
				shadow_grid[curr_i][curr_j][curr_k] += delta_s;
			}
		}
	}

	//printf("Shadows propogated.\n");

}

float TreeGenerator::get_light(int i, int j, int k)
{

	float s = shadow_grid[i][j][k];
	//printf("Shadow at (%d, %d, %d) = %f\n", i, j, k, s);
	float Q = c - s + a; //1.0 - 0.0+0.2
	//printf("Light at (%d, %d, %d) = %f\n", i, j, k, Q);
	if (Q < 0) Q = 0;
	return Q;

}

void TreeGenerator::min_voxel(int *min, std::vector<int*> voxels)
{
	//printf("Finding min voxel\n");
	//printf("\tvoxels.size() = %d\n", voxels.size());

	//shuffle
	std::random_shuffle(voxels.begin(), voxels.end());
	printf("first voxel = (%d, %d, %d)\n", voxels.at(0)[0], voxels.at(0)[1], voxels.at(0)[2]);

	float min_s = shadow_grid [voxels.at(0)[0]] [voxels.at(0)[1]] [voxels.at(0)[2]];
	min[0] = voxels.at(0)[0];
	min[1] = voxels.at(0)[1];
	min[2] = voxels.at(0)[2];

	int i;
	//printf("\tvoxels.size() = %d\n", voxels.size());
	for(i=0; i<voxels.size(); i++)
	{

		//printf("\t%d\n", i);
		float curr_s = shadow_grid [voxels.at(i)[0]] [voxels.at(i)[1]] [voxels.at(i)[2]];
		//printf("min_voxel:\tvoxel (%d, %d, %d) = %f\n", voxels.at(i)[0], voxels.at(i)[1], voxels.at(i)[2], curr_s);
		if(curr_s < min_s)
		{
			min_s = curr_s;
			min[0] = voxels.at(i)[0];
			min[1] = voxels.at(i)[1];
			min[2] = voxels.at(i)[2];
		}
	}
	//printf("min_voxel:\tmin voxel (%d, %d, %d) = %f\n", min[0], min[1], min[2], shadow_grid[min[0]][min[1]][min[2]]);
	//printf("Finished finding min voxel\n");

}

double TreeGenerator::vector_angle(vec3_t a, vec3_t b)
{

	//printf("vector_angle:\t a=<%f,%f,%f>\n", a[0], a[1], a[2]);
	//printf("vector_angle:\t b=<%f,%f,%f>\n", b[0], b[1], b[2]);

	double theta;
	double mag_a = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	double mag_b = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	//printf("vector_angle:\t mag_a = %f\n", mag_a);
	//printf("vector_angle:\t mag_b = %f\n", mag_b);
	double a_dot_b = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	//printf("vector_angle:\t a_dot_b = %f\n", a_dot_b);

	theta = acos(a_dot_b / (mag_a * mag_b));
	return theta;

}

//returns 15 points -> 3 layers of center and corner locations
void TreeGenerator::locations_from_voxel(int *voxel, vec3_t *locations)
{
	//y is easy
	locations[0][1] = locations[1][1] = locations[2][1] = locations[3][1] = locations[4][1] = voxel[1] + 0.1;
	locations[5][1] = locations[6][1] = locations[6][1] = locations[8][1] = locations[9][1] = voxel[1] + 0.5;
	locations[10][1] = locations[11][1] = locations[12][1] = locations[13][1] = locations[14][1] = voxel[1] + 0.9;

	//absolute center location (vox 0,0) is (99.5, 99.5)
	locations[0][0] = locations[5][0] = locations[10][0] = voxel[0] - (200-1)/2 + 0.5;
	locations[0][2] = locations[5][2] = locations[10][2] = voxel[2] - (200-1)/2 + 0.5;

	//rest can be found from this center location
	locations[1][0] = locations[6][0] = locations[11][0] = locations[0][0] + 0.4;
	locations[1][2] = locations[6][2] = locations[11][2] = locations[0][2] + 0.4;

	locations[2][0] = locations[7][0] = locations[12][0] = locations[0][0] + 0.4;
	locations[2][2] = locations[7][2] = locations[12][2] = locations[0][2] - 0.4;

	locations[3][0] = locations[8][0] = locations[13][0] = locations[0][0] - 0.4;
	locations[3][2] = locations[8][2] = locations[13][2] = locations[0][2] + 0.4;

	locations[4][0] = locations[9][0] = locations[14][0] = locations[0][0] - 0.4;
	locations[4][2] = locations[9][2] = locations[14][2] = locations[0][2] - 0.4;

}

void TreeGenerator::simple_location_from_voxel(int *voxel, vec3_t location)
{

	location[1] = voxel[1] + 0.5;
	location[0] = voxel[0] - (200-1)/2 + 0.5;
	location[2] = voxel[2] - (200-1)/2 + 0.5;

}

void TreeGenerator::normalize(vec3_t vec)
{
	double mag_vec = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	vec[0] /= mag_vec;
	vec[1] /= mag_vec;
	vec[2] /= mag_vec;
	//printf("normalize:\tvec = <%f, %f, %f>\n", vec[0], vec[1], vec[2]);
}

//for now, just assume perception_theta is always 90 degrees

void TreeGenerator::calculate_optimal_direction(node *start_node, double length, vec3_t dir)
{
	//printf("Determining direction\n");
	//printf("\t\tCurrent direction: (%f, %f, %f)\n", start_node->dir[0], start_node->dir[1], start_node->dir[2]);
	int start_voxel[3];
	voxel_from_location(start_voxel, start_node->location);

	std::vector<int*> neighbor_voxels;
	int center_voxel[3];
	vec3_t center_location = { start_node->location[0] + start_node->dir[0]*1.01,
							   start_node->location[1] + start_node->dir[1]*1.01,
							   start_node->location[2] + start_node->dir[2]*1.01 };
	voxel_from_location(center_voxel, center_location);
	//printf("Center voxel: [%d, %d, %d]\n", center_voxel[0], center_voxel[1], center_voxel[2]);

	//neighbor_voxels.push_back(center_voxel);

	//screw it. look at all neighboring voxels
	for(int i = -1; i<=1; i++)
	{
		for(int j = -1; j<=1; j++)
		{
			for(int k = -1; k<=1; k++)
			{
				if(i==0 && j==0 && k==0) 
					continue;
				int *curr_voxel = new int[3];
				curr_voxel[0] = start_voxel[0]+i;
				curr_voxel[1] = start_voxel[1]+j;
				curr_voxel[2] = start_voxel[2]+k;
				if(curr_voxel[0] >= 200 || curr_voxel[0] < 0 || curr_voxel[1] >= 200 || curr_voxel[1] < 0 || curr_voxel[2] >= 200 || curr_voxel[2] < 0)
					continue;
				
				vec3_t possible_locations[15];
				locations_from_voxel(curr_voxel, possible_locations);
				for(int m=0; m<15; m++)
				{
					vec3_t temp_dir = { possible_locations[m][0] - start_node->location[0],
										possible_locations[m][1] - start_node->location[1],
										possible_locations[m][2] - start_node->location[2] };
					double theta = vector_angle(temp_dir, start_node->dir);
					//printf("calculate_optimal_direction:\tvoxel(%d, %d, %d), theta = %f\n", curr_voxel[0], curr_voxel[1], curr_voxel[2], (theta*180/M_PI));
					if (theta > (45.0 * M_PI / 180.0))
					{
						//printf("calculate_optimal_direction:\tignoring\n");
					}
					else
					{
						//printf("calculate_optimal_direction:\tadding (%d, %d, %d)\n", curr_voxel[0], curr_voxel[1], curr_voxel[2]);
						neighbor_voxels.push_back(curr_voxel);
						break;
					}					
				}
				/*
				//printf("calculate_optimal_direction:\t voxel examining: [%d, %d, %d]\n", curr_voxel[0], curr_voxel[1], curr_voxel[2]);
				vec3_t possible_location;
				simple_location_from_voxel(curr_voxel, possible_location);
				vec3_t possible_dir = { possible_location[0] - start_node->location[0],
										possible_location[1] - start_node->location[1],
										possible_location[2] - start_node->location[2] };
				double theta = vector_angle(possible_dir, start_node->dir);
				if (theta > (45.0 * M_PI / 180.0))
					continue;
				else
				{
					neighbor_voxels.push_back(curr_voxel);
					break;
				}*/
			}
		}
	}

	//for(int i=0; i<neighbor_voxels.size(); i++)
//	{
		//printf("calculate_optimal_direction:\tneighbor_voxel[%d] = (%d, %d, %d)\n", i, neighbor_voxels.at(i)[0], neighbor_voxels.at(i)[1], neighbor_voxels.at(i)[2]);
//	}

	int min[3];

	min_voxel(min,neighbor_voxels);
	////printf("\t\tstart voxel: [%d, %d, %d]\n", start_voxel[0], start_voxel[1], start_voxel[2]);
	////printf("\t\tmin voxel: [%d, %d, %d]\n", min[0], min[1], min[2]);

	//after we find min voxel, figure out a possible location
	/*vec3_t end_location;
	vec3_t possible_locations[15];
	locations_from_voxel(min, possible_locations);
	for(int m=0; m<15; m++)
	{
		vec3_t temp_dir = { possible_locations[m][0] - start_node->location[0],
							possible_locations[m][1] - start_node->location[1],
							possible_locations[m][2] - start_node->location[2] };
		double theta = vector_angle(temp_dir, start_node->dir);
		if (theta > (45.0 * M_PI / 180.0))
			continue;
		else
		{
			end_location[0] = possible_locations[m][0];
			end_location[1] = possible_locations[m][1];
			end_location[2] = possible_locations[m][2];
			break;
		}					
	}*/

	vec3_t end_location;
	simple_location_from_voxel(min, end_location);

	dir[0] = end_location[0] - start_node->location[0];
	dir[1] = end_location[1] - start_node->location[1];
  	dir[2] = end_location[2] - start_node->location[2];

	////printf("\t\t dir = (%f, %f, %f)\n", dir[0], dir[1], dir[2]);
	normalize(dir);
	//tada
	////printf("\t\tdecided dir = (%f, %f, %f)\n", dir[0], dir[1], dir[2]);

}

void TreeGenerator::free_grid()
{

	int i, j;
	for(i=0; i<200; i++)
	{
		for(j=0; j<200; j++)
		{
			delete( shadow_grid[i][j] );
		}
	}

	delete( shadow_grid );

	//printf("Shadow grid freed\n");
}

void TreeGenerator::voxel_from_location(int *voxel, vec3_t location)
{
	//middle-bottom of voxel [200/2 - 1, 0, 200/2 - 1] is origin (location [0.0,0.0,0.0])
	//voxel width/height/length = 0.1

	//easy part is y
	voxel[1] = location[1];

	voxel[0] = location[0] + ((200-1)/2 + 0.5); 
	voxel[2] = location[2] + ((200-1)/2 + 0.5); 
}

void TreeGenerator::propogate_light_info()
{
	//start with most recent branches
	for(int i = branches.size() - 1; i >= 0; i--)
	{
		//printf("Branch %d\n....", i);
		double total_light = 0.0;
		int total_buds = 0;
		//aaaand start with the terminal node
		for(int j = branches.at(i).nodes.size()-1; j >= 0; j--)
		{
			//printf("\tNode %d\n....", j);
			int voxel[3];
			voxel_from_location(voxel, branches.at(i).nodes.at(j)->location);
			
			if(j > 0 && branches.at(i).nodes.at(j)->type != BRANCH)
			{
				//printf("\t\tNot a branch or base.\n");
				branches.at(i).nodes.at(j)->light_info = get_light(voxel[0], voxel[1], voxel[2]);
				branches.at(i).nodes.at(j)->avg_light_info = branches.at(i).nodes.at(j)->light_info;
				//printf("\t\tLight info = %f\n", branches.at(i).nodes.at(j)->light_info);
				//printf("\t\tAvg Light info = %f\n", branches.at(i).nodes.at(j)->avg_light_info);
				total_light += branches.at(i).nodes.at(j)->light_info;
				total_buds++;
			}
			else if(j > 0 && branches.at(i).nodes.at(j)->type == BRANCH)
			{
				//printf("\t\tIs a branch\n");
				branches.at(i).nodes.at(j)->light_info = branches.at(i).nodes.at(j)->branch_node->light_info;
				branches.at(i).nodes.at(j)->avg_light_info = branches.at(i).nodes.at(j)->branch_node->avg_light_info;
				//printf("\t\tLight info = %f\n", branches.at(i).nodes.at(j)->light_info);
				//printf("\t\tAvg Light info = %f\n", branches.at(i).nodes.at(j)->avg_light_info);
				total_light += branches.at(i).nodes.at(j)->light_info;
				total_buds += branches.at(i).nodes.at(j)->branch_node->total_buds;
			}
			else
			{
				//printf("\t\tIs a base.\n");
				branches.at(i).nodes.at(0)->light_info = total_light;
				//printf("\t\tTotal light = %f\n", total_light);
				branches.at(i).nodes.at(0)->total_buds = total_buds;
				branches.at(i).nodes.at(0)->avg_light_info = total_light / branches.at(i).nodes.at(0)->total_buds;
				//printf("\t\tAvg light info = %f\n", branches.at(i).nodes.at(0)->avg_light_info);
			}
		}
		//now sort the nodes on this branch (minus the base and the terminal nodes) by avg light info
		std::vector<node *> sorted;
		int num_nodes = branches.at(i).nodes.size()-2;
		//printf("Num nodes of branch %d: %d. Copying...\n", i, num_nodes);
		//copy references to nodes - ignore base and terminal
		for(int j=1; j<=num_nodes; j++)
		{
			//printf("Node %d\n", j);
			sorted.push_back(branches.at(i).nodes.at(j));
		}
		//good ole insertion sort
		//printf("Sorting...\n");
		for(int j=0; j<sorted.size(); j++)
		{
			node *key = sorted.at(j);
			int k;
			for(k = j-1; k >= 0; k--)
			{
				if(sorted.at(k)->avg_light_info > key->avg_light_info)
					break;
				sorted.at(k+1) = sorted.at(k);
			}
			sorted.at(k+1) = key;
		}
		//insert the terminal node at the beginning
		sorted.insert(sorted.begin(), branches.at(i).nodes.at( branches.at(i).nodes.size()-1 ));
		branches.at(i).sorted_nodes = sorted;
	}
}

double TreeGenerator::weight_function(int i, double w_min, double w_max, double kappa, int N)
{
	double weight = 0.0;
	if(i >= 0 && i < (kappa*N))
	{
		weight = - (w_min-w_max)/(kappa*N) * i + w_max;
	}
	else
		weight = w_min;


	return weight;
}

double TreeGenerator::get_resources(branch b, node *curr_node)
{

	//find sorted index of curr_node
	int sorted_indx = 5;
	for(int i=0; i<b.sorted_nodes.size(); i++)
	{
		if(b.sorted_nodes.at(i) == curr_node)
		{
			sorted_indx = i+1;
			break;
		}
	}
	//printf("sorted_index = %d\n", sorted_indx);

	//printf("\t\t\tGetting resources\n");
	//printf("Avg light = %f\n", curr_node->avg_light_info);
	double w_max = 1.0, w_min = 0.006, kappa = 0.5;

	double v_i = 0;
	double base_v = b.nodes.at(0)->base_v;
	//printf("\t\t\tbase_v = %f\n", base_v);
	double sum = 0.0;
	for(int j=1; j<=b.sorted_nodes.size(); j++)
	{
		double weight = weight_function(j, w_min, w_max, kappa, b.nodes.size()-1);
		//printf("\t\t\tSum ~ Avg light info = %f\n", b.nodes.at(j)->avg_light_info);
		//printf("\t\t\tSum ~ weight = %f\n", weight);
		sum += b.sorted_nodes.at(j-1)->avg_light_info * weight;
	}
	//printf("\t\t\tsum = %f\n", sum);
	double weight_i = weight_function(sorted_indx, w_min, w_max, kappa, b.nodes.size()-1);
	//printf("weight_i = %f\n", weight_i);
	v_i = base_v * (curr_node->avg_light_info * weight_i) / sum;
	//printf("v_i = %f\n", v_i);
	return v_i;
}

void TreeGenerator::calculate_location(node *start_node, node *end_node, double length)
{
	end_node->location[0] = start_node->location[0] + start_node->dir[0] * length;
	end_node->location[1] = start_node->location[1] + start_node->dir[1] * length;
	end_node->location[2] = start_node->location[2] + start_node->dir[2] * length;
}

/*
	node_type type;
	vec3_t location;
	vec3_t dir;
	int has_branch;
	node *branch_node;
	double light_info;
	double avg_light_info;
	int total_buds;
	double base_v;
	*/

void TreeGenerator::calculate_growth_direction(node *start_node, double length)
{
	vec3_t xi;							//optimal direction
	vec3_t eta = {0.0,1.0,0.0};			//tropism - grow upward toward light
	vec3_t def = { start_node->dir[0],	//default direction
				   start_node->dir[1],
				   start_node->dir[2] };
	calculate_optimal_direction(start_node, length, xi);

	//final direction is combination of three

	start_node->dir[0] = 0.5 * xi[0] + 0.3 * eta[0] + 1.0 * def[0];
	start_node->dir[1] = 0.5 * xi[1] + 0.3 * eta[1] + 1.0 * def[1];
	start_node->dir[2] = 0.5 * xi[2] + 0.3 * eta[2] + 1.0 * def[2];

	normalize(start_node->dir);

}


void TreeGenerator::grow_metamors(branch *b, node *start_node, int num_metamors, double length)
{
	//printf("Number of metamors in this step: %d\n", num_metamors);

	vec3_t start_dir;
	start_dir[0] = start_node->dir[0];
	start_dir[1] = start_node->dir[1];
	start_dir[2] = start_node->dir[2];

	for(int i=0; i<num_metamors; i++)
	{
		//node 0 was changed to axillary or base before this function was called
		if(i>0) start_node->type = AXILLARY_BUD;

		vec3_t start_dir;
		start_dir[0] = start_node->dir[0];
		start_dir[1] = start_node->dir[1];
		start_dir[2] = start_node->dir[2];
		int vox[3];
		voxel_from_location(vox, start_node->location);
		calculate_growth_direction(start_node, length);
		node *new_node = new node;
		new_node->type = TERMINAL_BUD;
		new_node->has_branch = false;
		new_node->dir[0] = start_node->dir[0];
		new_node->dir[1] = start_node->dir[1];
		new_node->dir[2] = start_node->dir[2];

		calculate_location(start_node, new_node, length);

		//add shadow values
		int curr_voxel[3];
		voxel_from_location(curr_voxel, new_node->location);
		propogate_shadow(curr_voxel[0], curr_voxel[1], curr_voxel[2]);

		//start_node is now the new node
		start_node = new_node;

		b->nodes.push_back(new_node);


	}

}


void TreeGenerator::get_bud_direction(vec3_t curr_dir, vec3_t bud_dir)
{
	//perhaps there is a better way to do this
	//but my brain can't math it

	//get start dir -> make y of curr_dir = 0
	vec3_t start_dir = {curr_dir[0], 0.0, curr_dir[2]};
	//printf("get_bud_direction:\t start_dir = <%f, %f, %f>\n", start_dir[0], start_dir[1], start_dir[2]);
	//if vector is about 0,0,0
	if(start_dir[0] < 0.02 && start_dir[1] < 0.02 && start_dir[2] < 0.02)
	{
		//then whatevs
		start_dir[0] = 1.0;
	}
	normalize(start_dir);
	//printf("get_bud_direction:\t start_dir = <%f, %f, %f>\n", start_dir[0], start_dir[1], start_dir[2]);

	//get angle between this start direction and the current direction
	double theta = vector_angle(start_dir, curr_dir);
	//printf("get_bud_direction:\t theta = %f\n", theta);
	//if curr_dir's y is less than 0, then adjust theta to be within 0 and 360deg
	if(curr_dir[1] < 0)
	{
		theta = 2*M_PI - theta;
	}

	//randomly decide to add -45 or 45deg
	if(rand()%2 == 0)
	{
		theta += M_PI / 4.0;
		//fix it if we went overboard
		if(theta > 2*M_PI)
			theta -= 2*M_PI;
	}
	else
	{
		theta -= M_PI / 4.0;
		if(theta < 0)
			theta = 2*M_PI + theta;
	}

	//then....
	//aaahhhh... ummmmmmm....
	//printf("get_bud_direction:\t new theta = %f\n", theta);
	bud_dir[1] = sin(theta);
	bud_dir[0] = start_dir[0] * cos(theta);
	bud_dir[2] = start_dir[2] * cos(theta);
	//printf("get_bud_direction:\t bud_dir = <%f, %f, %f>\n", bud_dir[0], bud_dir[1], bud_dir[2]);
	normalize(bud_dir);	//something like that?
	//printf("get_bud_direction:\t bud_dir = <%f, %f, %f>\n", bud_dir[0], bud_dir[1], bud_dir[2]);

}


void TreeGenerator::grow_tree_step()
{
	//printf("New step\n------------------------------------\n");
	//printf("Propogating light info...");
	propogate_light_info();
	//printf("Done!\n");
	double alpha = 2.0;
	double base_v = alpha * branches.at(0).nodes.at(0)->light_info;
	//printf("base_v = %f\n", base_v);
	branches.at(0).nodes.at(0)->base_v = base_v;

	//distribute resources, starting from base
	//printf("Distributing resources\n");
	int num_branches = branches.size();
	for(int i=0; i<num_branches; i++)
	{
		int num_nodes = branches.at(i).sorted_nodes.size();
		for(int j=0; j<num_nodes; j++)
		{
			node *curr_node = branches.at(i).sorted_nodes.at(j);
			double resources = get_resources(branches.at(i), curr_node);
			printf("RESOURCES = %f\n", resources);
			//printf("Branch %d, Node %d, resources = %f\n", i, j, resources);
			if(curr_node->type == BRANCH)
			{
				//printf("Setting branch base_v as %f\n", resources);
				curr_node->branch_node->base_v = resources;
			}
			else if(curr_node->type == AXILLARY_BUD)
			{
				int num_metamors = floor(resources);
				if(num_metamors == 0)
				{
					//printf("ignoring\n");
					continue;
				}

				//printf("creating a branch\n");
				double length = resources / num_metamors;
				curr_node->type = BRANCH;
				node *branch_node = new node;
				branch_node->location[0] = curr_node->location[0];
				branch_node->location[1] = curr_node->location[1];
				branch_node->location[2] = curr_node->location[2];

				//when creating a branch off a terminal bud, point (randomly) at +45deg or -45deg
				get_bud_direction(curr_node->dir, branch_node->dir);

				branch_node->base_v = resources;
				branch_node->type = BASE;
				curr_node->branch_node = branch_node;

				//create new branch
				branch *b = new branch;
				b->nodes.push_back(branch_node);
				grow_metamors(b, branch_node, num_metamors, length);
				branches.push_back(*b);
			}
			else if(curr_node->type == TERMINAL_BUD)
			{
				int num_metamors = floor(resources);
				if(num_metamors == 0)
				{
					//printf("ignoring\n");
					continue;
				}
				//printf("Expanding terminal bud\n");
				double length = resources / num_metamors;
				curr_node->type = AXILLARY_BUD;
				//use old branch
				grow_metamors(&branches.at(i), curr_node, num_metamors, length);
			}
		}
	}
}

void TreeGenerator::grow_tree(int num_steps)
{
	init_shadow_grid();
	init_tree();
	for(int i=0; i<num_steps; i++)
	{
		//printf("grow_tree: step %d\n", i);
		grow_tree_step();
	}
}

void TreeGenerator::print_tree()
{
	for(int i=0; i<branches.size(); i++)
	{
		//printf("\nTree branch #%d\n", i);
		//printf("------------------\n");
		for(int j=0; j<branches.at(i).nodes.size(); j++)
		{
			//printf("node %d location: (%f, %f, %f) - type: %d\n", j, branches.at(i).nodes.at(j)->location[0], branches.at(i).nodes.at(j)->location[1], branches.at(i).nodes.at(j)->location[2], branches.at(i).nodes.at(j)->type);
		}
	}
}


