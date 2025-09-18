import math
import numpy
import numpy.linalg


def normalizeColor(aColor):
    levels = 16
    factor = 255 / 16
    return (int(round(factor * round(aColor[0] * levels))),
            int(round(factor * round(aColor[1] * levels))),
            int(round(factor * round(aColor[2] * levels))))

class UnicolorMesh:
    def __init__(self):
        self.coord2index = {}
        self.faces = []

    def add(self, aVertices):
        indices = []
        for vertex in aVertices:
            if isinstance(vertex, numpy.ndarray):
                vertex = (vertex[0], vertex[1], vertex[2])
            if vertex in self.coord2index:
                index = self.coord2index[vertex]
            else:
                index = len(self.coord2index)
                self.coord2index[vertex] = index
            indices.append(index)
        self.faces.append(indices)

    def writeVertices(self, aFile, aColor):
        for vertex, _ in self.coord2index.items():
            aFile.write(f'{vertex[0]} {vertex[1]} {vertex[2]} {aColor[0]} {aColor[1]} {aColor[2]}\n')

    def writeFaces(self, aFile, aOffset):
        for face in self.faces:
            aFile.write(f'{len(face)} ')
            for index in face:
                aFile.write(f'{int(index + aOffset)} ')
            aFile.write('\n')

class StanfordPly:
    epsilon = 1e-7

    def __init__(self):
        self.meshes = {}

    # aColor already should be normalized.
    def add(self, aVertices, aColor):
        if not aColor in self.meshes:
            self.meshes[aColor] = UnicolorMesh()
        self.meshes[aColor].add(aVertices)

    def points3d(self, aX, aY, aZ, aColor, aRadius):
        color = normalizeColor(aColor)
        mmm = (aX - aRadius, aY - aRadius, aZ - aRadius)
        mmp = (aX - aRadius, aY - aRadius, aZ + aRadius)
        mpm = (aX - aRadius, aY + aRadius, aZ - aRadius)
        mpp = (aX - aRadius, aY + aRadius, aZ + aRadius)
        pmm = (aX + aRadius, aY - aRadius, aZ - aRadius)
        pmp = (aX + aRadius, aY - aRadius, aZ + aRadius)
        ppm = (aX + aRadius, aY + aRadius, aZ - aRadius)
        ppp = (aX + aRadius, aY + aRadius, aZ + aRadius)
        self.add((mmm, pmm, ppm, mpm), color)
        self.add((mmm, pmm, pmp, mmp), color)
        self.add((mpm, ppm, ppp, mpp), color)
        self.add((pmm, ppm, ppp, pmp), color)
        self.add((mmm, mpm, mpp, mmp), color)
        self.add((mmp, pmp, ppp, mpp), color)

    def plot3d(self, aX, aY, aZ, aColor, aRadius):
        nSegment = len(aX) - 1
        assert(nSegment >= 1)
        assert(nSegment == len(aY) - 1)
        assert(nSegment == len(aZ) - 1)
        color = normalizeColor(aColor)
        for i in range(nSegment):
            v0 = numpy.array([aX[i], aY[i], aZ[i]], numpy.float64)
            v1 = numpy.array([aX[i + 1], aY[i + 1], aZ[i + 1]], numpy.float64)
            displacement = v1 - v0
            length = numpy.linalg.norm(displacement)
            if length > StanfordPly.epsilon:
                displacement = displacement / length
                perpendicular0 = numpy.array([0, 0, 0], numpy.float64)
                if abs(displacement[1]) < StanfordPly.epsilon and abs(displacement[2]) < StanfordPly.epsilon:
                    perpendicular0[1] = 1.0
                    perpendicular0[2] = 0.0
                else:
                    denominator = math.sqrt(displacement[1] * displacement[1] + displacement[2] * displacement[2])
                    perpendicular0[1] = -displacement[2] / denominator
                    perpendicular0[2] =  displacement[1] / denominator
                perpendicular1 = numpy.cross(displacement, perpendicular0)
                mmm = v0 - perpendicular0 * aRadius
                mpp = v0 + perpendicular0 * aRadius
                mmp = v0 - perpendicular1 * aRadius
                mpm = v0 + perpendicular1 * aRadius
                pmm = v1 - perpendicular0 * aRadius
                ppp = v1 + perpendicular0 * aRadius
                pmp = v1 - perpendicular1 * aRadius
                ppm = v1 + perpendicular1 * aRadius
                self.add((mmm, pmm, ppm, mpm), color)
                self.add((mmm, pmm, pmp, mmp), color)
                self.add((mpm, ppm, ppp, mpp), color)
                self.add((pmm, ppm, ppp, pmp), color)
                self.add((mmm, mpm, mpp, mmp), color)
                self.add((mmp, pmp, ppp, mpp), color)

    def mesh(self, aX, aY, aZ, aColor):
        shape = aX.shape
        assert(shape == aY.shape and shape == aZ.shape)
        nI = shape[0] - 1
        nJ = shape[1] - 1
        assert(nI > 0 and nJ > 0)
        color = normalizeColor(aColor)
        for i in range(nI):
            for j in range(nJ):
                mm = numpy.array([aX[i, j], aY[i, j], aZ[i, j]], numpy.float64)
                mp = numpy.array([aX[i, j + 1], aY[i, j + 1], aZ[i, j + 1]], numpy.float64)
                pm = numpy.array([aX[i + 1, j], aY[i + 1, j], aZ[i + 1, j]], numpy.float64)
                pp = numpy.array([aX[i + 1, j + 1], aY[i + 1, j + 1], aZ[i + 1, j + 1]], numpy.float64)
                if numpy.linalg.norm(mm - mp) > StanfordPly.epsilon and numpy.linalg.norm(mp - pp) > StanfordPly.epsilon:
                    self.add((mm, mp, pp), color)
                if numpy.linalg.norm(pm - pp) > StanfordPly.epsilon and numpy.linalg.norm(pm - mm) > StanfordPly.epsilon:
                    self.add((pm, mm, pp), color)

    def triangularMesh(self, aX, aY, aZ, aTriangles, aColor):
        color = normalizeColor(aColor)
        assert(len(aX) == len(aY) and len(aY) == len(aZ))
        vertices = []
        for i in range(len(aX)):
            vertices.append((aX[i], aY[i], aZ[i]))
        for triangle in aTriangles:
            self.add((vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]]), color)

    def write(self, aName):
        nVertices = 0
        nFaces = 0
        for _, mesh in self.meshes.items():
            nVertices += len(mesh.coord2index)
            nFaces += len(mesh.faces)
        with open(aName, 'w', encoding='utf-8') as file:
            file.write('ply\nformat ascii 1.0\ncomment author: Balazs Bamer\ncomment object: RustBCA 3D visualization\n')
            file.write(f'element vertex {nVertices}\n')
            file.write('property float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\n')
            file.write(f'element face {nFaces}\n')
            file.write('property list uchar int vertex_indices\nend_header\n')

            for color, mesh in self.meshes.items():
                mesh.writeVertices(file, color)
            nVertices = 0
            for _, mesh in self.meshes.items():
                mesh.writeFaces(file, nVertices)
                nVertices += len(mesh.coord2index)
