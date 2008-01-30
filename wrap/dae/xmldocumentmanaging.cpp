#include <wrap/dae/xmldocumentmanaging.h>
#include <cassert>

XMLNode::XMLNode(XMLTag* tag)
:_tag(tag)
{
}

XMLNode::~XMLNode()
{
	delete _tag;
}

XMLLeafNode::XMLLeafNode(XMLLeafTag* leaftag)
:XMLNode(leaftag)
{
}

XMLLeafNode::~XMLLeafNode()
{

}


void XMLLeafNode::applyProcedure(Visitor& v)
{
	v(*this);
}

XMLInteriorNode::XMLInteriorNode(XMLTag* tag)
:XMLNode(tag)
{
}

XMLNode* XMLInteriorNode::son(int ii)
{
	assert((ii > 0) && (ii < _sons.size()));
	return _sons[ii];
}

QVector< XMLNode* > XMLInteriorNode::sons()
{
	return _sons;
}

XMLInteriorNode::~XMLInteriorNode()
{
	for(QVector< XMLNode* >::iterator it = _sons.begin();it != _sons.end();++it)
		delete (*it);
}

void XMLInteriorNode::applyProcedure(Visitor& v)
{
	v(*this);
}