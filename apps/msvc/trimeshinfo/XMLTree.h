
#include <utility>
#include <map>
#include <list>

#include <string>
#include <iostream>
#include "SlotsNode.h"
#include "ClassesNode.h"
#include "InstancesNode.h"


using namespace std;






class FacetNode: public Node
{
virtual void printNode();
virtual int qualifyNode();
};



class FacetsNode: public Node
{
	NodeGroup facets;
	virtual void printNode();
	virtual int qualifyNode();
};



class MainNode: public Node
{
public:

	MainNode(void){node_type = MAIN_NODE;};
	int node_type;
	list<pair<char*,char*> > headers;

	void addHeaders(char* str, char*val);
	virtual void printNode();
	virtual int qualifyNode();
};

void MainNode::addHeaders(char* str, char*val)
{
	headers.push_back(pair<char*,char*>(str,val));
}
void MainNode::printNode()
{

	cout<<"MainNode: node_type is "<<node_type<<"\n";
	list<pair<char*,char*> >::iterator it;
	for(it = headers.begin(); it!= headers.end(); ++it)
	{
		cout<<"MainNode: First element is "<< it->first<<"\n";
		cout<<"MainNode: Second element is "<<it->second<<"\n";
	}
}

int MainNode::qualifyNode()
{return node_type;}


class XMLTree
{
public:
	XMLTree(void){};
	~XMLTree(void){};
	NodeGroup root;
	

	// methods
	void initializeMain(char* xsn);
	void finalizeMain(char* xsn, char* val);
	void addHeaders(char* str, char*val);

	void addSlots(SlotNode* sn);
//	void addFacets();
	void addClasses(ClassNode* cn);
	void addInstances(InstanceNode* in);
	void printXMLTree();

};

void XMLTree::initializeMain(char* xsn)
{
	MainNode*	mn = new MainNode;
	
	mn->headers.push_back(pair<char*,char*>(xsn,""));
	root.Sons.push_back(mn);

}

void XMLTree::finalizeMain(char* xsn, char* val)
{
	MainNode*	mn = new MainNode;
	
	mn->headers.push_back(pair<char*,char*>(xsn,val));
	root.Sons.push_back(mn);

}


void XMLTree::addHeaders(char* str, char*val)
{
	MainNode* mn = (MainNode*) root.Sons.front();
	mn->headers.push_back(pair<char*,char*>(str,val));
}

void XMLTree::addSlots(SlotNode* sn)
{
	SlotsNode* sn0 = new SlotsNode; // 1 solo

	sn0->addSlot(sn);
	root.Sons.push_back(sn0);
	

}

void XMLTree::addClasses(ClassNode* cn)
{
	ClassesNode* cn0 = new ClassesNode; // 1 solo

	cn0->addClass(cn);
	root.Sons.push_back(cn0);
}

void XMLTree::addInstances(InstanceNode* in)
{
	InstancesNode* in0 = new InstancesNode; // 1 solo

	in0->addInstance(in);
	root.Sons.push_back(in0);
}

void XMLTree::printXMLTree()
{	
	string ext,s("XMLFile");
	cout<<"Preparing to create XML file"<<"\n";
	cout<<"enter the name of the file "<< endl;
	cin >> ext;


	s.append(ext);
	s.append(".xml");

	const char* filename = s.c_str();
	FILE* fp = fopen(filename, "w");

	list<Node*>::iterator it;
	list<Node*>::iterator it2;
	list<Node*>::iterator it3;
	list<pair<char*,char*> >::iterator lit;
	MainNode* mn;
	SlotsNode* sns;
	SlotNode* sn;
	OwnSlotNode* osn;	
	ClassesNode* csn;	
	ClassNode* cn;	
	InstancesNode* isn;	
	InstanceNode* in;	
	int nn = 0;
	for(it = root.Sons.begin(); it!=root.Sons.end(); ++it){
		cout<<"Creating Node #"<< nn<<"\n";
		cout<<"Node Type is "<< (*it)->qualifyNode()<<"\n";
		switch((*it)->qualifyNode())
		{
		case MAIN_NODE:
			mn = (MainNode*)(*it);
			fprintf(fp,"<");
			for(lit = mn->headers.begin(); lit!= mn->headers.end(); ++lit)
					fprintf(fp,"%s%s", lit->first,lit->second );
				fprintf(fp,"> \n");
			
			break;
		
		case SLOTS_NODE:
			sns = (SlotsNode*)(*it);
			fprintf(fp,"<slots> \n");

			for(it2 = sns->slot.Sons.begin(); it2!=sns->slot.Sons.end(); ++it2)
			{
				sn = (SlotNode*) (*it2);
				fprintf(fp,"\t<slot>");
				for(it3 = sn->own_slot.Sons.begin(); it3!=sn->own_slot.Sons.end(); ++it3)
				{
					osn = (OwnSlotNode*) (*it3);
					fprintf(fp,"<own-slot slot-name=\":%s\">\n",osn->name);
					fprintf(fp,"\t\t<entry type=\"%s\">\n",osn->entry.type);
					fprintf(fp,"\t\t\t<value>\n");
					fprintf(fp,"\t\t\t%s\n",osn->entry.value.value);
					fprintf(fp,"\t\t\t</value>\n");
					fprintf(fp,"\t\t</entry>\n");
					fprintf(fp,"\t</own-slot>");
				}
				fprintf(fp,"</slot>\n");
			}
			fprintf(fp,"</slots>\n");

			break;
	case CLASSES_NODE:
			csn = (ClassesNode*)(*it);
			fprintf(fp,"<classes> \n");

			for(it2 = csn->classn.Sons.begin(); it2!=csn->classn.Sons.end(); ++it2)
			{
				cn = (ClassNode*) (*it2);
				fprintf(fp,"\t<class>");
				fprintf(fp,"<own-slots>\n");
				for(it3 = cn->own_slots.own_slot.Sons.begin(); it3!=cn->own_slots.own_slot.Sons.end(); ++it3)
				{
					osn = (OwnSlotNode*) (*it3);
					fprintf(fp,"\t\t<own-slot slot-name=\":%s\">\n",osn->name);
					fprintf(fp,"\t\t\t<entry type=\"%s\">\n",osn->entry.type);
					fprintf(fp,"\t\t\t\t<value>\n");
					fprintf(fp,"\t\t\t%s\n",osn->entry.value.value);
					fprintf(fp,"\t\t\t\t</value>\n");
					fprintf(fp,"\t\t\t</entry>\n");
					fprintf(fp,"\t\t</own-slot>\n");
				}
				fprintf(fp,"\t</own-slots>\n");
				fprintf(fp,"</class>\n");
			}
			fprintf(fp,"</classes>\n");

			break;
	case INSTANCES_NODE:
			isn = (InstancesNode*)(*it);
			fprintf(fp,"<instances> \n");

			for(it2 = isn->instances.Sons.begin(); it2!=isn->instances.Sons.end(); ++it2)
			{
				in = (InstanceNode*) (*it2);
				fprintf(fp,"\t<instance>\n");
				fprintf(fp,"\t\t<id>\n");
				fprintf(fp,"\t\t%s\n", in->id);
				fprintf(fp,"\t\t</id>\n");
				fprintf(fp,"\t\t<type>\n");
				fprintf(fp,"\t\t%s\n", in->type);
				fprintf(fp,"\t\t</type>\n");
				fprintf(fp,"\t\t<own-slots>\n");
				for(it3 = in->own_slots.own_slot.Sons.begin(); it3!=in->own_slots.own_slot.Sons.end(); ++it3)
				{
					osn = (OwnSlotNode*) (*it3);
					fprintf(fp,"\t\t\t<own-slot slot-name=\":%s\">\n",osn->name);
					fprintf(fp,"\t\t\t\t<entry type=\"%s\">\n",osn->entry.type);
					fprintf(fp,"\t\t\t\t<value>\n");
					fprintf(fp,"\t\t\t\t%s\n",osn->entry.value.value);
					fprintf(fp,"\t\t\t\t</value>\n");
					fprintf(fp,"\t\t\t\t</entry>\n");
					fprintf(fp,"\t\t\t</own-slot>\n");
				}
				fprintf(fp,"\t\t</own-slots>\n");
				fprintf(fp,"\t</instance>\n");
			}
			fprintf(fp,"</instances>\n");

			break;
		}
		++nn;
	}
	fclose(fp);
}