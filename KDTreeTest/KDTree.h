#ifndef __KDTREE_H__
#define __KDTREE_H__


#define MAX_DIM  20             //max space dimension
#define KD_MAX_POINTS 2000000    //max KD tree points number  
#define SD 10					//current space dimension
#define RMAX 1000
#define POINTS_NUM 10000
#define TEST_NUM 50

#define KDTEMPLATE	template <class Xtype>
#define KDNODE		KDNode<Xtype>
#define KDTREE		KDTree<Xtype>

KDTEMPLATE inline Xtype distance2(Xtype* x1, Xtype* x2, int dim)
{
	Xtype S = 0;
	for(int k=0; k < dim; k++)
		S+= (x1[k]-x2[k])*(x1[k]-x2[k]);
	return S;
}


KDTEMPLATE	bool equal(Xtype* x1, Xtype* x2, int dim)
{
	for(int k=0; k < dim; k++)
	{
		if(x1[k] != x2[k])
			return false;
	}

	return true ;
}

//KDTreeNode template class implementation
KDTEMPLATE class KDNode
{
//member functions
public:
	int axis ;
	Xtype x[SD];
	unsigned int id ;
	bool checked ;
	bool orientation ;

	KDNode(Xtype* x0, int axis0);

	KDNODE*	Insert(Xtype* x);
	KDNODE*	FindParent(Xtype* x0);

	KDNODE* Parent ;
	KDNODE* Left ;
	KDNODE* Right ;
};


KDTEMPLATE
KDNODE::KDNode(Xtype* x0, int axis0)
{
	axis = axis0;
	for(int k=0; k<SD; k++)
		x[k] = x0[k];

	Left = Right = Parent = NULL ;
	checked = false ;
	id = 0;
}


KDTEMPLATE
KDNODE*	KDNODE::FindParent(Xtype* x0)
{
	KDNODE* parent ;
	KDNODE* next = this ;
	int split ;
	while(next)
	{
		split = next->axis  ;
		parent = next ;
		if(x0[split] > next->x[split])
			next = next->Right ;
		else
			next = next->Left ;
	}
	return parent ;
}

KDTEMPLATE
KDNODE*	KDNODE::Insert(Xtype* p)
{
	KDNODE* parent = FindParent(p);
	if(equal(p, parent->x, SD))
		return NULL ;

	KDNODE* newNode = new KDNODE(p, parent->axis +1 < SD? parent->axis+1:0);
	newNode->Parent = parent ;

	if(p[parent->axis] > parent->x[parent->axis])
	{
		parent->Right = newNode ;
		newNode->orientation = 1 ;
	}
	else
	{
		parent->Left = newNode ;
		newNode->orientation = 0 ;
	}

	return newNode ;
}



KDTEMPLATE
class KDTree
{
public:
	KDNODE*  Root ;
	KDTree();

	bool				add(Xtype* x);
	KDNODE*				find_nearest(Xtype* x0);
	KDNODE*				find_nearest_brute(Xtype* x) ;

	inline	void		check_subtree(KDNODE* node, Xtype* x);
	inline  void		set_bounding_cube(KDNODE* node, Xtype* x);
	inline KDNODE*		search_parent(KDNODE* parent, Xtype* x);
	void				uncheck();

	LARGE_INTEGER		TimeStart, TimeFinish; 
	LARGE_INTEGER		CounterFreq;  

	void GetTimerFrequency(){	
		QueryPerformanceFrequency(&CounterFreq);
	}

	void StartTimer(){
	QueryPerformanceCounter(&TimeStart);
	}

	void StopTimer(){
	QueryPerformanceCounter(&TimeFinish);
	}

	double GetElapsedTime(){
		double calc_time = (double)(TimeFinish.LowPart - TimeStart.LowPart)/(double)CounterFreq.LowPart;
		return 1000*calc_time ;
	}

public:	
	Xtype				d_min ;
	KDNODE*				nearest_neighbour ;

	int					KD_id  ;

	KDNODE*				List[KD_MAX_POINTS] ;
	int					nList ;

	KDNODE*				CheckedNodes[KD_MAX_POINTS] ;
	int					checked_nodes ;

	Xtype				x_min[SD], x_max[SD]; 
	bool				max_boundary[SD], min_boundary[SD];
	int					n_boundary ;

};

KDTEMPLATE
KDTREE::KDTree()
{
	Root = NULL ;
	KD_id = 1;
	nList = 0;
}


KDTEMPLATE
bool KDTree<Xtype>::add(Xtype* x)
{
	if(nList >= KD_MAX_POINTS-1)
		return 0; //can't add more points
	
	if(!Root)
	{
		Root =  new KDNODE(x,0);
		Root->id = KD_id++ ;
		List[nList++] = Root ;
	}
	else
	{
		KDNODE* pNode ;
		if((pNode=Root->Insert(x)))
		{
			pNode->id = KD_id++ ;
			List[nList++] = pNode;
		}
	}

	return true ;
}


//sequential nearest neighbour search
KDTEMPLATE	
KDNODE* KDTree<Xtype>::find_nearest_brute(Xtype* x) 
{
	if(!Root)
		return NULL;

	KDNODE* nearest = Root ;
	Xtype d ;
	d_min = distance2(Root->x, x, SD);
	for(int n=1; n<nList; n++)
	{
		d =  distance2(List[n]->x, x, SD);
		if(d < d_min)
		{
			nearest = List[n] ;
			d_min = d;
		}
	}

	return nearest ;
}


KDTEMPLATE
KDNODE* KDTREE::find_nearest(Xtype* x)
{
	if(!Root)
		return NULL ;

	checked_nodes = 0;

	KDNODE* parent = Root->FindParent(x);
	nearest_neighbour = parent ;
	d_min = distance2(x, parent->x, SD); ;

	if(equal(x, parent->x, SD))
		return nearest_neighbour ;

	search_parent(parent, x);
	uncheck();

	return nearest_neighbour ;
}




KDTEMPLATE
void KDTREE::check_subtree(KDNODE* node, Xtype* x)
{
	if(!node || node->checked)
		return ;

	CheckedNodes[checked_nodes++] = node ;
	node->checked = true ;
	set_bounding_cube(node, x);
	
	int dim = node->axis ;
	Xtype d= node->x[dim] - x[dim];

	if (d*d > d_min) {
		if (node->x[dim] > x[dim])
		  check_subtree(node->Left, x);
		else 
		  check_subtree(node->Right, x);
	}
	// If the distance from the key to the current value is 
	// less than the nearest distance, we still need to look
	// in both directions.
	else {
		check_subtree(node->Left, x);
		check_subtree(node->Right, x);
	}
}

KDTEMPLATE
void KDTREE::set_bounding_cube(KDNODE* node, Xtype* x)
{
	if(!node)
		return ;
	int d=0;
	Xtype dx;
	for(int k=0; k<SD; k++)
	{
		dx = node->x[k]-x[k];
		if(dx > 0)
		{
			dx *= dx ;
			if(!max_boundary[k])
			{
				if(dx > x_max[k])
					x_max[k] = dx;
				if(x_max[k]>d_min)
				{
					max_boundary[k] =true;
					n_boundary++ ;
				}
			}
		}
		else 
		{
			dx *= dx ;
			if(!min_boundary[k])
			{
				if(dx > x_min[k])
					x_min[k] = dx;
				if(x_min[k]>d_min)
				{
					min_boundary[k] =true;
					n_boundary++ ;
				}
			}
		}
		d+=dx;
		if(d>d_min)
			return ;

	}
	
	if(d<d_min)
	{
		d_min = d;
		nearest_neighbour = node ;
	}
}

KDTEMPLATE
KDNODE* KDTREE::search_parent(KDNODE* parent, Xtype* x)
{
	for(int k=0; k<SD; k++)
	{
		x_min[k] = x_max[k] = 0;
		max_boundary[k] = min_boundary[k] = 0;
	}
	n_boundary = 0;

	Xtype dx;
	KDNODE* search_root = parent ;
	while(parent && n_boundary != 2*SD)
	{	
		check_subtree(parent, x);
		search_root = parent ;
		parent = parent->Parent ;
	}

	return search_root ;
}

KDTEMPLATE
void KDTREE::uncheck()
{
	for(int n=0; n<checked_nodes; n++)
		CheckedNodes[n]->checked = false;
}

#endif