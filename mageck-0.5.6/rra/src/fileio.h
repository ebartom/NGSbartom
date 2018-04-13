#ifndef FILEIO_H
#define FILEIO_H

#include "classdef.h"

//Read input file. File Format: <item id> <group id> <list id> <value>. Return 1 if success, -1 if failure
int ReadFile(char *fileName, GROUP_STRUCT*& groups, int maxGroupNum, int *groupNum, LIST_STRUCT *lists, int maxListNum, int *listNum);

//Save group information to output file. Format <group id> <number of items in the group> <lo-value> <false discovery rate>
int SaveGroupInfo(char *fileName, GROUP_STRUCT *groups, int groupNum);





#endif
