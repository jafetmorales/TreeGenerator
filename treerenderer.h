#ifndef TREERENDERER_H
#define TREERENDERER_H

#include "treegenerator.h"

typedef struct branch_model
{
	int branchid;
	GLuint vbo;
	float *verts;
	int num_verts;
} branch_model;

typedef struct tree_model
{
	branch_model **branches;
	int size;
} tree_model;

branch_model *get_model(branch b);
tree_model *get_tree_model(std::vector<branch>);
void init_tree(tree_model *tm);
void render_tree(tree_model *tm);

#endif