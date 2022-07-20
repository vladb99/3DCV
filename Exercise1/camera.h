#ifndef CAMERA_H
#define CAMERA_H

#include <QVector3D>
#include <QMatrix4x4>
#include <cube.h>

class Camera
{
public:
    Camera(QVector3D _position, QVector3D _rotation, float _f);
    ~Camera();

    QMatrix4x4 getRotationMatrix();
    QVector3D getImagePrinciplePoint();
    QMatrix4x4 getTransformationMatrix();
    QVector3D affine_projection(QVector3D point_to_project);
    QVector3D homogeneous_projection(QVector3D point_to_project);
    QVector3D explicit_projection(QVector3D point_to_project);
    std::vector<QVector3D> project(Cube cube, int algorithm);
    QVector3D move_into_plane(QVector2D projected_point);

private:
    QVector3D position;
    QVector3D rotation;
    float f;

public:
    QVector3D getPosition() const { return position; }
    QVector3D getRotation() const { return rotation; }
    float getF() const { return f; }
};

#endif // CAMERA_H
