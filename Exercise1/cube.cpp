#include "cube.h"

Cube::Cube(QVector3D _position, QVector3D _rotation, float _size) :
    position(_position),
    rotation(_rotation),
    size(_size)
{}

Cube::~Cube()
{}

Cube::Cube(QVector3D _position, QVector3D _position_2){
    position = _position;
    position_2 = _position_2;
}

std::vector<QVector3D> Cube::get_cube_from_2_points(){
    std::vector<QVector3D> tmp;

    QVector4D a = QVector4D(position_2.x()/ 0.05f, position_2.y()/ 0.05f, position.z()/ 0.05f,       1);
    QVector4D b = QVector4D(position.x()/ 0.05f,   position_2.y()/ 0.05f, position.z()/ 0.05f,       1);
    QVector4D c = QVector4D(position.x()/ 0.05f,   position.y()/ 0.05f,   position.z()/ 0.05f,       1);
    QVector4D d = QVector4D(position_2.x()/ 0.05f, position.y()/ 0.05f,   position.z()/ 0.05f,       1);
    QVector4D e = QVector4D(position_2.x()/ 0.05f, position_2.y()/ 0.05f, position_2.z()/ 0.05f,     1);
    QVector4D f = QVector4D(position.x()/ 0.05f,   position_2.y()/ 0.05f, position_2.z()/ 0.05f,     1);
    QVector4D g = QVector4D(position.x()/ 0.05f,   position.y()/ 0.05f,   position_2.z()/ 0.05f,     1);
    QVector4D h = QVector4D(position_2.x()/ 0.05f, position.y()/ 0.05f,   position_2.z()/ 0.05f,     1);

    //A-B
    tmp.push_back(a.toVector3D());
    tmp.push_back(b.toVector3D());

    //A-D
    tmp.push_back(a.toVector3D());
    tmp.push_back(d.toVector3D());

    //C-B
    tmp.push_back(c.toVector3D());
    tmp.push_back(b.toVector3D());

    //C-D
    tmp.push_back(c.toVector3D());
    tmp.push_back(d.toVector3D());

    //A-E
    tmp.push_back(a.toVector3D());
    tmp.push_back(e.toVector3D());

    //B-F
    tmp.push_back(b.toVector3D());
    tmp.push_back(f.toVector3D());

    //C-G
    tmp.push_back(c.toVector3D());
    tmp.push_back(g.toVector3D());

    //D-H
    tmp.push_back(d.toVector3D());
    tmp.push_back(h.toVector3D());

    //E-F
    tmp.push_back(e.toVector3D());
    tmp.push_back(f.toVector3D());

    //E-H
    tmp.push_back(e.toVector3D());
    tmp.push_back(h.toVector3D());

    //G-F
    tmp.push_back(g.toVector3D());
    tmp.push_back(f.toVector3D());

    //G-H
    tmp.push_back(g.toVector3D());
    tmp.push_back(h.toVector3D());

    return tmp;
}

std::vector<QVector3D> Cube::get_cube_points(){
     std::vector<QVector3D> tmp;
     QVector4D translate = QVector4D(position.x(), position.y(), position.z(), 0);

     float angleX =  rotation.x() * M_PI / 180;
     float angleY =  rotation.y() * M_PI / 180;
     float angleZ =  rotation.z() * M_PI / 180;

     QVector4D a = QVector4D(-size, size, 0, 1);
     QVector4D b = QVector4D(0, size, 0, 1);
     QVector4D c = QVector4D(0, 0, 0, 1);
     QVector4D d = QVector4D(-size, 0, 0, 1);
     QVector4D e = QVector4D(-size, size, -size, 1);
     QVector4D f = QVector4D(0, size, -size, 1);
     QVector4D g = QVector4D(0, 0, -size, 1);
     QVector4D h = QVector4D(-size, 0, -size, 1);

     // Rotation
     a = rotate_point_x_axis(a, angleX);
     b = rotate_point_x_axis(b, angleX);
     c = rotate_point_x_axis(c, angleX);
     d = rotate_point_x_axis(d, angleX);
     e = rotate_point_x_axis(e, angleX);
     f = rotate_point_x_axis(f, angleX);
     g = rotate_point_x_axis(g, angleX);
     h = rotate_point_x_axis(h, angleX);

     a = rotate_point_y_axis(a, angleY);
     b = rotate_point_y_axis(b, angleY);
     c = rotate_point_y_axis(c, angleY);
     d = rotate_point_y_axis(d, angleY);
     e = rotate_point_y_axis(e, angleY);
     f = rotate_point_y_axis(f, angleY);
     g = rotate_point_y_axis(g, angleY);
     h = rotate_point_y_axis(h, angleY);

     a = rotate_point_z_axis(a, angleZ);
     b = rotate_point_z_axis(b, angleZ);
     c = rotate_point_z_axis(c, angleZ);
     d = rotate_point_z_axis(d, angleZ);
     e = rotate_point_z_axis(e, angleZ);
     f = rotate_point_z_axis(f, angleZ);
     g = rotate_point_z_axis(g, angleZ);
     h = rotate_point_z_axis(h, angleZ);

     //translation

     a =  a + translate;
     b =  b + translate;
     c =  c + translate;
     d =  d + translate;
     e =  e + translate;
     f =  f + translate;
     g =  g + translate;
     h =  h + translate;

    //A-B
    tmp.push_back(a.toVector3D());
    tmp.push_back(b.toVector3D());

    //A-D
    tmp.push_back(a.toVector3D());
    tmp.push_back(d.toVector3D());

    //C-B
    tmp.push_back(c.toVector3D());
    tmp.push_back(b.toVector3D());

    //C-D
    tmp.push_back(c.toVector3D());
    tmp.push_back(d.toVector3D());

    //A-E
    tmp.push_back(a.toVector3D());
    tmp.push_back(e.toVector3D());

    //B-F
    tmp.push_back(b.toVector3D());
    tmp.push_back(f.toVector3D());

    //C-G
    tmp.push_back(c.toVector3D());
    tmp.push_back(g.toVector3D());

    //D-H
    tmp.push_back(d.toVector3D());
    tmp.push_back(h.toVector3D());

    //E-F
    tmp.push_back(e.toVector3D());
    tmp.push_back(f.toVector3D());

    //E-H
    tmp.push_back(e.toVector3D());
    tmp.push_back(h.toVector3D());

    //G-F
    tmp.push_back(g.toVector3D());
    tmp.push_back(f.toVector3D());

    //G-H
    tmp.push_back(g.toVector3D());
    tmp.push_back(h.toVector3D());

    return tmp;
}

QVector4D Cube::rotate_point_y_axis(QVector4D vector, float angle) {
    QMatrix4x4 rotation_matrix = QMatrix4x4(cos(angle), 0, sin(angle), 0, 0, 1, 0, 0, -sin(angle), 0, cos(angle), 0, 0,0,0,1 );
    return rotation_matrix * vector;
}

QVector4D Cube::rotate_point_x_axis(QVector4D vector, float angle) {
    QMatrix4x4 rotation_matrix = QMatrix4x4( 1, 0, 0, 0, 0, cos(angle), -sin(angle), 0, 0, sin(angle), cos(angle), 0, 0,0,0,1);
    return rotation_matrix * vector;
}

QVector4D Cube::rotate_point_z_axis(QVector4D vector, float angle) {
    QMatrix4x4 rotation_matrix = QMatrix4x4(cos(angle), -sin(angle), 0, 0, sin(angle), cos(angle), 0, 0, 0, 0, 1, 0,0,0,0,1);
    return rotation_matrix * vector;
}
