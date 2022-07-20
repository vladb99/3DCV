#ifndef CUBE_H
#define CUBE_H


#include <QVector3D>
#include <QMatrix4x4>
#include <QColor>

class Cube
{
public:
    QVector3D position;
    QVector3D position_2;
    Cube(QVector3D _position, QVector3D _rotation, float _size);
    Cube(QVector3D _position, QVector3D _position_2);
    ~Cube();

    std::vector<QVector3D> get_cube_points();    
    std::vector<QVector3D> get_cube_from_2_points();

    QVector4D rotate_point_x_axis(QVector4D vector, float angle);
    QVector4D rotate_point_y_axis(QVector4D vector, float angle);
    QVector4D rotate_point_z_axis(QVector4D vector, float angle);

private:

    QVector3D rotation;
    float size;

public:
    QVector3D getPosition() const { return position; }
    QVector3D getRotation() const { return rotation; }
    float getSize() const { return size; }
};
#endif // CUBE_H
