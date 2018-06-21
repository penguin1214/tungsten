#include "MeshIO.hpp"

#include <stdio.h>

#include "FileUtils.hpp"
#include "ObjLoader.hpp"
#include "IntTypes.hpp"

#include <tinyformat/tinyformat.hpp>

namespace Tungsten {

namespace MeshIO {

bool saveObj(const Path &path, const std::vector<Vertex> &verts, const std::vector<TriangleI> &tris);

bool loadWo3(const Path &path, std::vector<Vertex> &verts, std::vector<TriangleI> &tris)
{

    InputStreamHandle stream = FileUtils::openInputStream(path);
    if (!stream)
        return false;

    uint64 numVerts, numTris;
    FileUtils::streamRead(stream, numVerts);
    verts.resize(size_t(numVerts));
    FileUtils::streamRead(stream, verts);
    FileUtils::streamRead(stream, numTris);
    tris.resize(size_t(numTris));
    FileUtils::streamRead(stream, tris);

#if 0
	Path p("wo3output.obj");
	saveObj(p, verts, tris);
#endif

    return true;
}

bool saveWo3(const Path &path, const std::vector<Vertex> &verts, const std::vector<TriangleI> &tris)
{
    OutputStreamHandle stream = FileUtils::openOutputStream(path);
    if (!stream)
        return false;

    FileUtils::streamWrite(stream, uint64(verts.size()));
    FileUtils::streamWrite(stream, verts);
    FileUtils::streamWrite(stream, uint64(tris.size()));
    FileUtils::streamWrite(stream, tris);

    return true;
}

std::vector<Vec3f> loadGradient(std::string fname) {
	std::cout << "loading gradient... " << std::endl;
	const char *fn = fname.c_str();
	FILE *fp = fopen(fn, "rb");
	if (fp == NULL) {
		perror("Error: ");
	}
	long bufsize;
	if (fseek(fp, 0L, SEEK_END) == 0) {
		if (bufsize == -1) { perror("Error: "); }
		bufsize = ftell(fp);
		// read back to file start
		if (fseek(fp, 0L, SEEK_SET) != 0) { perror("Error: "); }
	}

	int s = bufsize / (3 * sizeof(float));
	std::vector<Vec3f> gradient;
	// gradient.reserve(bufsize);
	//fread(&gradient[0], 3*sizeof(float), bufsize, fp);
	float x, y, z;
	for (int i = 0; i < s; i++) {
		fread(&x, sizeof(float), 1, fp);
		fread(&y, sizeof(float), 1, fp);
		fread(&z, sizeof(float), 1, fp);
		gradient.push_back(Vec3f(x, y, z).normalized());
	}
	return gradient;
}

std::vector<Vec3f> loadVec3Bin(std::string fname, int gridx, int gridy, int gridz) {
	const char *fn = fname.c_str();
	FILE *fp = fopen(fn, "rb");
	if (fp == NULL) {
		perror("Error");
	}

	fread(&gridx, sizeof(int), 1, fp);
	fread(&gridy, sizeof(int), 1, fp);
	fread(&gridz, sizeof(int), 1, fp);

	std::vector<Vec3f> gradient;
	gradient.reserve(gridx*gridy*gridz);

	int idx;
	float x, y, z;
	for (int i = 0; i < gridz; i++) {
		for (int j = 0; j < gridy; j++) {
			for (int k = 0; k < gridx; k++) {
				idx = i*gridx*gridy + j*gridx + k;
				fread(&x, sizeof(float), 1, fp);
				fread(&y, sizeof(float), 1, fp);
				fread(&z, sizeof(float), 1, fp);
				gradient[i] = Vec3f(x, y, z);
			}
		}
	}
	fclose(fp);
	std::cout << "gradient loaded." << std::endl;
	return gradient;
}

bool loadObj(const Path path, std::vector<Vertex> &verts, std::vector<TriangleI> &tris)
{
	bool ret = ObjLoader::loadGeometryOnly(path, verts, tris);
	
#if 0
	Path p("test.wo3");
	saveWo3(p, verts, tris);
#endif
	return ret;
}

bool saveObj(const Path &path, const std::vector<Vertex> &verts, const std::vector<TriangleI> &tris)
{
    OutputStreamHandle stream = FileUtils::openOutputStream(path);
    if (!stream)
        return false;

    for (const Vertex &v : verts)
        tfm::format(*stream, "v %f %f %f\n", v.pos().x(), v.pos().y(), v.pos().z());
    for (const Vertex &v : verts)
        tfm::format(*stream, "vn %f %f %f\n", v.normal().x(), v.normal().y(), v.normal().z());
    for (const Vertex &v : verts)
        tfm::format(*stream, "vt %f %f\n", v.uv().x(), v.uv().y());
    for (const TriangleI &t : tris)
        tfm::format(*stream, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
            t.v0 + 1, t.v0 + 1, t.v0 + 1,
            t.v1 + 1, t.v1 + 1, t.v1 + 1,
            t.v2 + 1, t.v2 + 1, t.v2 + 1);

    return true;
}

bool load(const Path &path, std::vector<Vertex> &verts, std::vector<TriangleI> &tris)
{
    if (path.testExtension("wo3"))
        return loadWo3(path, verts, tris);
    else if (path.testExtension("obj"))
        return loadObj(path, verts, tris);
    return false;
}

bool save(const Path &path, const std::vector<Vertex> &verts, const std::vector<TriangleI> &tris)
{
    if (path.testExtension("wo3"))
        return saveWo3(path, verts, tris);
    else if (path.testExtension("obj"))
        return saveObj(path, verts, tris);
    return false;
}

}

}
