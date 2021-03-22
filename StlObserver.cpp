//
// Created by Nikita Kruk on 18.03.21.
//

#include "StlObserver.hpp"

#include <string>
#include <fstream>
#include <array>

StlObserver::StlObserver()
{

}

StlObserver::~StlObserver()
{

}

void StlObserver::SaveCellMesh(const CellMesh &cell_mesh) const
{
  std::string file_name("/Users/nikita/Documents/Projects/DeformableCellModel/Cell_Mesh.stl");
  std::ofstream file(file_name, std::ios::out | std::ios::trunc | std::ios::binary);

  // UINT8[80]    – Header                 -     80 bytes
  std::vector<char> header(80, '_');
  file.write((char *) &header[0], header.size() * sizeof(char));

  // UINT32       – Number of triangles    -      4 bytes
  unsigned int number_of_triangles = cell_mesh.GetNumFaces();
  file.write((char *) &number_of_triangles, sizeof(unsigned int));

  // foreach triangle                      - 50 bytes:
  //     REAL32[3] – Normal vector             - 12 bytes
  //     REAL32[3] – Vertex 1                  - 12 bytes
  //     REAL32[3] – Vertex 2                  - 12 bytes
  //     REAL32[3] – Vertex 3                  - 12 bytes
  //     UINT16    – Attribute byte count      -  2 bytes
  // end
  auto faces = cell_mesh.GetFaces();
  const auto &face_normals = cell_mesh.GetNormalsForFaces();
  auto nodes = cell_mesh.GetNodes();
  std::array<float, kDim> face_normal_to_save{0.0f}, vertex_to_save{0.0f};
  for (int f = 0; f < faces.size(); ++f)
  {
    face_normal_to_save[0] = face_normals[f][0];
    face_normal_to_save[1] = face_normals[f][1];
    face_normal_to_save[2] = face_normals[f][2];
    file.write((char *) &face_normal_to_save[0], kDim * sizeof(float));

    for (int vertex_idx : faces[f])
    {
      vertex_to_save[0] = nodes[vertex_idx][0];
      vertex_to_save[1] = nodes[vertex_idx][1];
      vertex_to_save[2] = nodes[vertex_idx][2];
      file.write((char *) &vertex_to_save[0], kDim * sizeof(float));
    } // vertex_idx

    static const unsigned short attribute_byte_count = 0;
    file.write((char *) &attribute_byte_count, sizeof(attribute_byte_count));
  } // f

  file.close();
}