#include <OpenGL/gl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "treerenderer.h"



//for now setting arbitrary radius of 0.2, number of points 5

struct branch_model *get_model(branch b)
{
	int num_points = 5;
	float radius = 0.2;

	float bottom[num_points][3];
	float top[num_points][3];

	float sin_vals[num_points];
	float cos_vals[num_points];

	branch_model *b_model = (struct branch_model *) malloc(sizeof(struct branch_model));
	b_model->verts = (float *) malloc( (b.nodes.size()) * (num_points*2 + 2) * 3 * sizeof(float) );
	//b_model->branchid = b->branchid;

	int i;

	//if i learned anything in comp org/arch it's that sin and cos are expensive
	//since theyre the same for every point, let's not recalculate them each time
	for(i=0; i<num_points; i++)
	{
		sin_vals[i] = sin( (i+1) * (2*M_PI/num_points) );
		cos_vals[i] = cos( (i+1) * (2*M_PI/num_points) );

		//well while we're here we might as well get the first bottom set of points
		bottom[i][0] = b.nodes.at(0)->location[0] + cos_vals[i]*radius;
		bottom[i][1] = b.nodes.at(0)->location[1] + sin_vals[i]*radius;
		bottom[i][2] = b.nodes.at(0)->location[2];
	}

	//printf("num allocated: %d\n", (b.nodes.size()) * (num_points*2 + 2));
	//printf("b->size-1=%d\n", b.nodes.size());
	//printf("num_points*2+2=%d\n", num_points*2+2);


	int vert_indx = 0;
	for(i=0; i<b.nodes.size()-1; i++)
	{
		////printf("------BOTTOM------\n");
		int j;
		for(j=0; j<num_points; j++)
		{
		//	//printf("(%f, %f, %f)\n", bottom[j][0], bottom[j][1], bottom[j][2]);
		}
		////printf("-------TOP-------\n");
		for(j=0; j < (num_points+1); j++)
		{
			if(j < num_points)
			{
				top[j][0] = b.nodes.at(i+1)->location[0] + cos_vals[j]*radius;
				top[j][1] = b.nodes.at(i+1)->location[1] + sin_vals[j]*radius;
				top[j][2] = b.nodes.at(i+1)->location[2];
				////printf("(%f, %f, %f)\n", top[j][0], top[j][1], top[j][2]);
			}

			//b_model->verts[vert_indx] = (float *) malloc(3*sizeof(float));
			//memcpy(b_model->verts[vert_indx], bottom[j%num_points], 3*sizeof(float));
			b_model->verts[vert_indx++] = bottom[j%num_points][0];
			b_model->verts[vert_indx++] = bottom[j%num_points][1];
			b_model->verts[vert_indx++] = bottom[j%num_points][2];

			////printf("vertex %d: (%f, %f, %f)\n", vert_indx, b_model->verts[vert_indx][0], b_model->verts[vert_indx][1], b_model->verts[vert_indx][2]);
			////printf("vert #%d\n", vert_indx+1);
			//vert_indx++;

			//b_model->verts[vert_indx] = (float *) malloc(3*sizeof(float));
			//memcpy(b_model->verts[vert_indx], top[j%num_points], 3*sizeof(float));
			b_model->verts[vert_indx++] = top[j%num_points][0];
			b_model->verts[vert_indx++] = top[j%num_points][1];
			b_model->verts[vert_indx++] = top[j%num_points][2];
			////printf("vert #%d\n", vert_indx+1);
			////printf("vertex %d: (%f, %f, %f)\n", vert_indx, b_model->verts[vert_indx][0], b_model->verts[vert_indx][1], b_model->verts[vert_indx][2]);
			//vert_indx++;
			//
			if(j > 0)
				memcpy(bottom[j%num_points], top[j%num_points], 3*sizeof(float));
		}

	}

	b_model->num_verts = (vert_indx-1)/3;
	//printf("vert_indx-1=%d\n", vert_indx-1);
/*
	for(i=0; i<vert_indx-1; i+=3)
	{
		//printf("vertex %d: (%f, %f, %f)\n", i, b_model->verts[i], b_model->verts[i+1], b_model->verts[i+2]);
	}
*/


	//glGenBuffers(1, &b_model->vbo);
	//glBindBuffer(GL_ARRAY_BUFFER, b_model->vbo);
	//glBufferData(GL_ARRAY_BUFFER, b_model->num_verts*3*sizeof(float), b_model->verts, GL_STATIC_DRAW);

	//printf("done\n");

	return b_model;
}

struct tree_model *get_tree_model(std::vector<branch> branches)
{
	tree_model *tm;
	tm = (tree_model *) malloc(sizeof(tree_model));
	tm->size = branches.size();
	tm->branches = (branch_model **) malloc(tm->size * sizeof(branch_model *));
	
	int i;
	for(i=0; i<tm->size; i++)
	{
		tm->branches[i] = get_model(branches.at(i));
	}
	return tm;
}



void init_branch(branch_model *bm)
{
	glGenBuffers(1, &bm->vbo);
	glBindBuffer(GL_ARRAY_BUFFER, bm->vbo);
	glBufferData(GL_ARRAY_BUFFER, bm->num_verts*3*sizeof(float), bm->verts, GL_STATIC_DRAW);

}

void init_tree(tree_model *tm)
{
	int i=0;
	for(i=0; i<tm->size; i++)
	{
		init_branch(tm->branches[i]);
	}

}

void render_branch(branch_model *bm, int i)
{

	glColor3f(0.58,0.29,0.0);
	glBindBuffer(GL_ARRAY_BUFFER, bm->vbo);
	glVertexPointer(3, GL_FLOAT, 0, NULL);
	glEnableClientState(GL_VERTEX_ARRAY);

	glDrawArrays(GL_TRIANGLE_STRIP, 0, bm->num_verts);
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void render_tree(tree_model *tm)
{
	int i;
	for(i=0; i<tm->size; i++)
	{
		render_branch(tm->branches[i],i);
	}
}


