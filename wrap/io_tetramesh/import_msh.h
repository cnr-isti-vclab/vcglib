#ifndef __VCGLIB_IMPORTTETMSH_H
#define __VCGLIB_IMPORTTETMSH_H

#include <iostream>

namespace vcg
{
namespace tetra
{
namespace io
{
template <class MeshType>
class ImporterMSH
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::TetraType TetraType;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::TetraIterator TetraIterator;

    enum ErrorCodes
    {
        INVALID_FORMAT = 1,
        INVALID_VERSION,
        NOT_IMPLEMENTED,
        IO_ERROR
    };

    static inline void parseWhiteSpace(std::ifstream &fin)
    {
        //we don't want to consume non whitespace bytes, just peek it..
        char next = fin.peek();

        while (next == '\n' || next == ' ' || next == '\t' || next == '\r')
        {
            fin.get();
            next = fin.peek();
        }
    }

    static int parseNodes(MeshType &m, std::ifstream &fin, bool binary)
    {
        int numOfNodes;
        fin >> numOfNodes;
        if (numOfNodes < 0)
            return INVALID_FORMAT;

        VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertices(m, numOfNodes);

        if (binary)
        {
            size_t lineBytes = (4 + 3 * 8); //int index + 3 * double coords
            size_t bytes = numOfNodes * lineBytes;
            char *data = new char[bytes];
            parseWhiteSpace(fin);
            fin.read(data, bytes);

            for (int i = 0; i < numOfNodes; ++i)
            {
                int index = *reinterpret_cast<int *>(&data[i * lineBytes]) - 1;
                if (index < 0)
                    return INVALID_FORMAT;

                m.vert[index].P().X() = *reinterpret_cast<float *>(&data[i * lineBytes + 4]);
                m.vert[index].P().Y() = *reinterpret_cast<float *>(&data[i * lineBytes + 4 + 8]);
                m.vert[index].P().Z() = *reinterpret_cast<float *>(&data[i * lineBytes + 4 + 2 * 8]);
            }
            delete[] data;
        }
        else
        {
            for (int i = 0; i < numOfNodes; ++i)
            {
                int index;
                fin >> index;
                --index;

                if (index < 0)
                    return INVALID_FORMAT;

                fin >> m.vert[index].P().X();
                fin >> m.vert[index].P().Y();
                fin >> m.vert[index].P().Z();
            }
        }
        return 0;
    }

    static int parseElements(MeshType &m, std::ifstream &fin, bool binary)
    {
        int numOfElements;
        fin >> numOfElements;
        if (numOfElements < 0)
            return INVALID_FORMAT;

        TetraIterator ti = vcg::tri::Allocator<MeshType>::AddTetras(m, numOfElements);

        if (binary)
        {
            parseWhiteSpace(fin);
            size_t parsedElems = 0;

            while (parsedElems < numOfElements)
            {
                int index, type, tags;
                fin.read((char *)&index, sizeof(int));
                fin.read((char *)&type, sizeof(int));
                fin.read((char *)&tags, sizeof(int));
                --index;

                //check for tetra type
                if (type != 4)
                    return NOT_IMPLEMENTED;
                //check index validity
                if (index < 0)
                    return INVALID_FORMAT;

                //read tags and throw them
                for (size_t j = 0; j < tags; ++j)
                {
                    int tag;
                    fin.read((char *)&tag, sizeof(int));
                }

                //read element nodes
                TetraType t = m.tetra[index];
                for (int i = 0; i < 4; ++i)
                {
                    int nodeIndex;
                    fin.read((char *)&nodeIndex, sizeof(int));
                    --nodeIndex;

                    if (nodeIndex < 0 || nodeIndex >= m.VN())
                        return INVALID_FORMAT;

                    t.V(i) = &m.vert[nodeIndex];
                }
                ++parsedElems;
            }
        }
        else
        {
            for (int i = 0; i < numOfElements; ++i)
            {
                int index, type, tags;
                fin >> index >> type >> tags;
                --index;

                // std::cerr << index << std::endl;

                //check for tetra type
                if (type != 4)
                    return NOT_IMPLEMENTED;
                //check index validity
                if (index < 0)
                    return INVALID_FORMAT;
                //read tags and throw them
                for (size_t j = 0; j < tags; ++j)
                {
                    int tag;
                    fin >> tag;
                }

                TetraType * t = &m.tetra[index];
                for (int i = 0; i < 4; ++i)
                {
                    int nodeIndex;
                    fin >> nodeIndex;
                    --nodeIndex;

                    if (nodeIndex < 0 || nodeIndex > m.VN())
                        return INVALID_FORMAT;
                    // std::cerr << nodeIndex << std::endl;

                    t->V(i) = &m.vert[nodeIndex];
                }
            }
        }

        return 0;
    }

    static int parseMshMesh(MeshType &m, std::string &filename)
    {
        std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);

        if (!fin.is_open())
            return IO_ERROR;

        std::string lookAhead;

        fin >> lookAhead;
        if (lookAhead != "$MeshFormat")
            return INVALID_FORMAT;


        double version;
        int type, dataSize;

        fin >> version >> type >> dataSize;

        if (version != 2.2)
            return INVALID_VERSION;

        bool binary = (type == 1);

        if (dataSize != 8)
            return INVALID_FORMAT;

        // Read in extra info from binary header.
        if (binary)
        {
            int one;
            parseWhiteSpace(fin);
            fin.read(reinterpret_cast<char *>(&one), sizeof(int));
            if (one != 1)
            {
                std::cerr << "Warning: binary msh file " << filename
                          << " is saved with different endianness than this machine."
                          << std::endl;
                throw NOT_IMPLEMENTED;
            }
        }

        lookAhead.clear();
        fin >> lookAhead;
        if (lookAhead != "$EndMeshFormat")
            return INVALID_FORMAT;

        std::cerr << "reading nodes and elements" << std::endl;

        while (!fin.eof())
        {
            lookAhead.clear();
            fin >> lookAhead;

            if (lookAhead == "$Nodes")
            {
                int res = parseNodes(m, fin, binary);

                std::cerr << "reading nodes: " << res << std::endl;


                if (res != 0)
                    return res;
                std::cerr << "reading nodes" << std::endl;

                fin >> lookAhead;
                if (lookAhead != "$EndNodes")
                    return INVALID_FORMAT;
            }
            else if (lookAhead == "$Elements")
            {
                int res = parseElements(m, fin, binary);

                if (res != 0)
                    return res;

                std::cerr << "elements" << std::endl;


                fin >> lookAhead;
                if (lookAhead != "$EndElements")
                    return INVALID_FORMAT;
            }
            else if (lookAhead == "$NodeData")
            {
                return NOT_IMPLEMENTED;
                // parse_node_field(fin);
                // fin >> lookAhead;
                // if (lookAhead != "$EndNodeData")
                //     return INVALID_FORMAT;
            }
            else if (lookAhead == "$ElementData")
            {
                return NOT_IMPLEMENTED;
                // parse_element_field(fin);
                // fin >> lookAhead;
                // if (lookAhead != "$EndElementData")
                //     return INVALID_FORMAT;
            }
            else if (fin.eof())
            {
                break;
            }
            else
            {
                return INVALID_FORMAT;
                // parse_unknown_field(fin, lookAhead);
            }
        }
        fin.close();
        return 0;
    }

  public:
    static int Open(MeshType &m, const char *filename, CallBackPos *cb = 0)
    {
        std::string name(filename);
        return parseMshMesh(m, name);
    }
};
} // namespace io
} // namespace tetra
} // namespace vcg

#endif
