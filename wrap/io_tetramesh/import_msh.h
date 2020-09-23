#ifndef __VCGLIB_IMPORTTETMSH_H
#define __VCGLIB_IMPORTTETMSH_H

#include <iostream>

namespace vcg {
namespace tetra {
namespace io {

template <class MeshType>
class MshInfo
{
    typedef std::map<std::string, double *> FieldMap;

    FieldMap nodeFields;
    FieldMap elemFields;


    template <typename Scalar, int Dimensions>
    struct AttribTraits
    {
        typedef std::vector<Scalar> Type;
        enum { Dimension = Dimensions };
    };

    template <typename Scalar>
    struct AttribTraits<Scalar, 1>
    {
            typedef Scalar Type;
            enum { Dimension = 1 };
    };

    template <typename AttrHandle, typename Scalar, int Dimension>
    struct AttribHelper
    {
        static void assign(AttrHandle & handle, int index, const double * data)
        {
            for (int i=0; i<Dimension; ++i)
            {
                handle[index][i] = Scalar(data[Dimension * index + i]);
            }
        }
    };

    template <typename AttrHandle, typename Scalar>
    struct AttribHelper<AttrHandle, Scalar, 1>
    {
        static void assign(AttrHandle & handle, int index, const double * data)
        {
            handle[index] = Scalar(data[index]);
        }
    };

    template <bool PerNode, typename Scalar>
    static void fillMeshAttributes(const std::string & attrib_name, int attrib_dim, MeshType & mesh, const double * data)
    {
        if (PerNode)
        {
            switch(attrib_dim)
            {
            case 1 : fillMeshWithAttributePerNode<Scalar, 1>(attrib_name, mesh, data); break;
            case 2 : fillMeshWithAttributePerNode<Scalar, 2>(attrib_name, mesh, data); break;
            case 3 : fillMeshWithAttributePerNode<Scalar, 3>(attrib_name, mesh, data); break;
            case 4 : fillMeshWithAttributePerNode<Scalar, 4>(attrib_name, mesh, data); break;
            case 5 : fillMeshWithAttributePerNode<Scalar, 5>(attrib_name, mesh, data); break;
            case 6 : fillMeshWithAttributePerNode<Scalar, 6>(attrib_name, mesh, data); break;
            case 7 : fillMeshWithAttributePerNode<Scalar, 7>(attrib_name, mesh, data); break;
            case 8 : fillMeshWithAttributePerNode<Scalar, 8>(attrib_name, mesh, data); break;
            case 9 : fillMeshWithAttributePerNode<Scalar, 9>(attrib_name, mesh, data); break;
            default : throw std::string("Dimension of custom attribute vector unsupported");
            }
        }
        else
        {
            switch(attrib_dim)
            {
            case 1 : fillMeshWithAttributePerElement<Scalar, 1>(attrib_name, mesh, data); break;
            case 2 : fillMeshWithAttributePerElement<Scalar, 2>(attrib_name, mesh, data); break;
            case 3 : fillMeshWithAttributePerElement<Scalar, 3>(attrib_name, mesh, data); break;
            case 4 : fillMeshWithAttributePerElement<Scalar, 4>(attrib_name, mesh, data); break;
            case 5 : fillMeshWithAttributePerElement<Scalar, 5>(attrib_name, mesh, data); break;
            case 6 : fillMeshWithAttributePerElement<Scalar, 6>(attrib_name, mesh, data); break;
            case 7 : fillMeshWithAttributePerElement<Scalar, 7>(attrib_name, mesh, data); break;
            case 8 : fillMeshWithAttributePerElement<Scalar, 8>(attrib_name, mesh, data); break;
            case 9 : fillMeshWithAttributePerElement<Scalar, 9>(attrib_name, mesh, data); break;
            default : throw std::string("Dimension of custom attribute vector unsupported");
            }
        }
    }

    template <typename Scalar, int Dimension/* = 1*/>
    static void fillMeshWithAttributePerNode(const std::string & attrib_name, MeshType & mesh, const double * data)
    {
        typedef typename AttribTraits<Scalar, Dimension>::Type AttrType;
        typedef typename MeshType::template PerVertexAttributeHandle<AttrType> AttrHandle;

        AttrHandle handle = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<AttrType>(mesh, attrib_name);
        size_t num_nodes = size_t(mesh.VN());

        for (int i=0; i<int(num_nodes); ++i)
            AttribHelper<AttrHandle, Scalar, Dimension>::assign(handle, i, data);
    }

    template <typename Scalar, int Dimension/* = 1*/>
    static void fillMeshWithAttributePerElement(const std::string & attrib_name, MeshType & mesh, const double * data)
    {
        typedef typename AttribTraits<Scalar, Dimension>::Type AttrType;
        typedef typename MeshType::template PerFaceAttributeHandle<AttrType> AttrHandle;

        AttrHandle handle = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<AttrType>(mesh, attrib_name);
        size_t num_elements = size_t(mesh.TN());

        for (int i=0; i<int(num_elements); ++i)
            AttribHelper<AttrHandle, Scalar, Dimension>::assign(handle, i, data);
    }
};

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

                m.vert[index].P().X() = *reinterpret_cast<double *>(&data[i * lineBytes + 4]);
                m.vert[index].P().Y() = *reinterpret_cast<double *>(&data[i * lineBytes + 4 + 8]);
                m.vert[index].P().Z() = *reinterpret_cast<double *>(&data[i * lineBytes + 4 + 2 * 8]);
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
                //MSH in binary format has a elem-header 3*4 bytes: {elems_type, numElems, tagsPerElem}
                //followed by the list of elems under this header and eventually a new header and list.
                int type, elements, tags;
                fin.read((char *)&type,     sizeof(int));
                fin.read((char *)&elements, sizeof(int));
                fin.read((char *)&tags,     sizeof(int));

                //check for tetra type
                if (type != 4)
                    return NOT_IMPLEMENTED;

                //read tags and throw them
                for (size_t j = 0; j < tags; ++j)
                {
                    int tag;
                    fin.read((char *)&tag, sizeof(int));
                }

                //foreach element
                for (int i = 0; i < elements; ++i)
                {
                    int index;
                    fin.read((char *)&index, sizeof(int));
                    --index;

                    //check index validity
                    if (index < 0)
                        return INVALID_FORMAT;

                    //read element nodes
                    TetraType * t = &m.tetra[index];
                    for (int i = 0; i < 4; ++i)
                    {
                        int nodeIndex;
                        fin.read((char *)&nodeIndex, sizeof(int));
                        --nodeIndex;

                        if (nodeIndex < 0 || nodeIndex >= m.VN())
                            return INVALID_FORMAT;

                        t->V(i) = &m.vert[nodeIndex];
                    }
                    ++parsedElems;
                }
            }
        }
        else
        {
            for (int i = 0; i < numOfElements; ++i)
            {
                int index, type, tags;
                fin >> index >> type >> tags;
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
                    fin >> tag;
                }

                TetraType *t = &m.tetra[index];
                for (int i = 0; i < 4; ++i)
                {
                    int nodeIndex;
                    fin >> nodeIndex;
                    --nodeIndex;

                    if (nodeIndex < 0 || nodeIndex > m.VN())
                        return INVALID_FORMAT;

                    t->V(i) = &m.vert[nodeIndex];
                }
            }
        }

        return 0;
    }

    static int parseDataField(MeshType &m, std::ifstream &fin, bool binary)
    {
        int numString, numReal, numInteger;

        fin >> numString;

        std::string *strTags = new std::string[numString];
        for (int i = 0; i < numString; ++i)
        {
            parseWhiteSpace(fin);
            fin >> strTags[i];
        }

        fin >> numReal;

        double *doubleTags = new double[numReal];
        for (int i = 0; i < numReal; ++i)
            fin >> doubleTags[i];

        fin >> numInteger;

        if (numString <= 0 || numInteger < 3)
            return INVALID_FORMAT;

        int *integerTags = new int[numInteger];

        for (int i = 0; i < numInteger; ++i)
            fin >> integerTags[i];

        std::string fieldName = strTags[0];
        int fieldComponents = integerTags[1];
        int fieldSize = integerTags[2];

        double *fieldVec = new double[fieldComponents * fieldSize];

        delete[] strTags;
        delete[] doubleTags;
        delete[] integerTags;

        if (binary)
        {
            size_t totalBytes = (4 + 8 * fieldComponents) * fieldSize;
            char *data = new char[totalBytes];
            parseWhiteSpace(fin);
            fin.read(data, totalBytes);

            for (int i = 0; i < fieldSize; ++i)
            {
                int index = *reinterpret_cast<int *>(&data[i * (4 + fieldComponents * 8)]);
                --index;

                if (index < 0)
                    return INVALID_FORMAT;

                //values
                int baseIndex = i * (4 + fieldComponents * 8) + 4;

                for (int j = 0; j < fieldComponents; ++j)
                    fieldVec[index * fieldComponents + j] = *reinterpret_cast<float *>(&data[baseIndex + j * 8]);
            }
        }
        else
        {
            for (int i = 0; i < fieldSize; ++i)
            {
                int index;
                fin >> index;
                --index;

                if (index < 0)
                    return INVALID_FORMAT;

                for (int j = 0; j < fieldComponents; ++j)
                    fin >> fieldVec[index * fieldComponents + j];
            }
        }
    }

    static int parseNodeData(MeshType &m, MshInfo<MeshType> & info, std::ifstream &fin, bool binary)
    {
        return parseDataField(m, fin, binary);
    }

    static int parseElementData(MeshType &m, MshInfo<MeshType> & info, std::ifstream &fin, bool binary)
    {
        return parseDataField(m, fin, binary);
    }

    static int parseUnsupportedTag(std::ifstream &fin, std::string &tag)
    {
        std::cerr << "found unsupported tag" << std::endl;
        std::string tagName = tag.substr(1, tag.size() - 1);

        std::string tagEnd = tag.substr(0, 1) + "End" + tagName;

        std::string buf;
        while (buf != tagEnd && !fin.eof())
            fin >> buf;
        return 0;
    }

    static int parseMshMesh(MeshType &m, std::string &filename, MshInfo<MeshType> & info)
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

        // Read endiannes info in binary header...it's a 1 used to detect endiannes.
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

        while (!fin.eof())
        {
            lookAhead.clear();
            fin >> lookAhead;

            if (lookAhead == "$Nodes")
            {
                int res = parseNodes(m, fin, binary);
                if (res != 0)
                    return res;

                fin >> lookAhead;
                if (lookAhead != "$EndNodes")
                    return INVALID_FORMAT;
            }
            else if (lookAhead == "$Elements")
            {
                int res = parseElements(m, fin, binary);

                if (res != 0)
                    return res;

                fin >> lookAhead;
                if (lookAhead != "$EndElements")
                    return INVALID_FORMAT;
            }
            else if (lookAhead == "$NodeData")
            {
                parseNodeData(m, info, fin, binary);
                fin >> lookAhead;
                if (lookAhead != "$EndNodeData")
                    return INVALID_FORMAT;
            }
            else if (lookAhead == "$ElementData")
            {
                parseElementData(m, info, fin, binary);
                fin >> lookAhead;
                if (lookAhead != "$EndElementData")
                    return INVALID_FORMAT;

            }
            else if (fin.eof())
            {
                break;
            }
            else
            {
                parseUnsupportedTag(fin, lookAhead);
            }
        }
        fin.close();
        return 0;
    }

public:
    static int Open(MeshType &m, const char *filename, CallBackPos *cb = 0)
    {
        MshInfo<MeshType> info;
        return Open(m, filename, info, cb);
    }

    static int Open(MeshType &m, const char *filename, MshInfo<MeshType> & info, CallBackPos *cb = 0)
    {
        std::string name(filename);
        return parseMshMesh(m, name, info);
    }
};
} // namespace io
} // namespace tetra
} // namespace vcg

#endif
