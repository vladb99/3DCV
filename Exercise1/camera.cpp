#include "camera.h"

Camera::Camera(QVector3D _position, QVector3D _rotation, float _f) :
    position(_position),
    rotation(_rotation),
    f(_f)
{}

Camera::~Camera()
{}

QMatrix4x4 Camera::getRotationMatrix() {
    QMatrix4x4 rotation_matrix_x;
    QMatrix4x4 rotation_matrix_y;
    QMatrix4x4 rotation_matrix_z;

    float x =  rotation.x() * M_PI / 180;
    float y =  rotation.y() * M_PI / 180;
    float z =  rotation.z() * M_PI / 180;

    rotation_matrix_x = QMatrix4x4( 1, 0, 0, 0, 0, cos(x), -sin(x), 0, 0, sin(x), cos(x), 0, 0,0,0,1);
    rotation_matrix_y = QMatrix4x4(cos(y), 0, sin(y), 0, 0, 1, 0, 0, -sin(y), 0, cos(y), 0, 0,0,0,1 );
    rotation_matrix_z = QMatrix4x4(cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0,0,0,0,1);

    return rotation_matrix_z * rotation_matrix_y * rotation_matrix_x;
}

QVector3D Camera::getImagePrinciplePoint() {
    return (getTransformationMatrix() * QVector4D(0, 0, f, 1)).toVector3D();
}

QMatrix4x4 Camera::getTransformationMatrix() {
    QMatrix4x4 rotation_and_translation_matrix = getRotationMatrix();
    rotation_and_translation_matrix.setColumn(3, QVector4D(position.x(), position.y(), position.z(), 1));
    return rotation_and_translation_matrix;
}

QVector3D Camera::move_into_plane(QVector2D projected_point) {
    QVector3D H = getImagePrinciplePoint();
    QVector2D diff = projected_point - QVector2D(H.x(), H.y());
    QVector3D direction = (getRotationMatrix() * QVector4D(diff.x(), diff.y(), 0, 1)).toVector3D();
    return H + direction;
}

// Works!
QVector3D Camera::affine_projection(QVector3D point_to_project) {
    QVector3D N = position;
    QVector3D H = getImagePrinciplePoint();
    QVector3D P = point_to_project;

    float focal_length = f;

    QMatrix4x4 camera_rotation_matrix_t = getRotationMatrix().transposed();
    QVector3D diff_PN = P - N;
    QVector4D rotated_point = camera_rotation_matrix_t * QVector4D(diff_PN.x(), diff_PN.y(), diff_PN.z(), 1);
    QVector2D projected_point_without_z = QVector2D(H.x(), H.y()) - (-focal_length / rotated_point.z()) * QVector2D(rotated_point.x(), rotated_point.y());

    return move_into_plane(projected_point_without_z);
}

// Works!
QVector3D Camera::homogeneous_projection(QVector3D point_to_project) {
    QVector3D N = position;
    QVector3D H = getImagePrinciplePoint();
    QVector3D P = point_to_project;

    float focal_length = f;

    QMatrix4x4 camera_rotation_matrix_t = getRotationMatrix().transposed();
    QVector4D rotated_N = -camera_rotation_matrix_t * QVector4D(N.x(), N.y(), N.z(), 1);
    QMatrix4x4 middle_matrix = camera_rotation_matrix_t;
    middle_matrix.setColumn(3, rotated_N);
    QVector4D result1 = middle_matrix * QVector4D(P.x(), P.y(), P.z(), 1);
    QMatrix4x4 last_matrix = QMatrix4x4(focal_length, 0, H.x(), 0,
                                        0, focal_length, H.y(), 0,
                                        0, 0, 1, 0, // -1 to 1 in this row fixes the problem!
                                        0, 0, 0, 1);
    QVector3D result2 = (last_matrix * result1).toVector3D();
    result2 = result2 / result2.z();
    QVector2D projected_point_without_z = QVector2D(result2.x(), result2.y());

    return move_into_plane(projected_point_without_z);
}

QVector3D Camera::explicit_projection(QVector3D point_to_project) {
    // Should work
    QMatrix4x4 camera_rotation_matrix_t = getRotationMatrix().transposed();
    float r11 = camera_rotation_matrix_t(0, 0);
    float r12 = camera_rotation_matrix_t(1, 0);
    float r13 = camera_rotation_matrix_t(2, 0);
    float r21 = camera_rotation_matrix_t(0, 1);
    float r22 = camera_rotation_matrix_t(1, 1);
    float r23 = camera_rotation_matrix_t(2, 1);
    float r31 = camera_rotation_matrix_t(0, 2);
    float r32 = camera_rotation_matrix_t(1, 2);
    float r33 = camera_rotation_matrix_t(2, 2);

    QVector3D N = position;
    QVector3D H = getImagePrinciplePoint();
    QVector3D P = point_to_project;

    float focal_length = f;

    float x_row_0 = r11*(P.x()-N.x()) + r21*(P.y() - N.y()) + r31*(P.z() - N.z());
    float x_row_1 = r13*(P.x()-N.x()) + r23*(P.y() - N.y()) + r33*(P.z() - N.z());

    float y_row_0 = r12*(P.x()-N.x()) + r22*(P.y() - N.y()) + r32*(P.z() - N.z());
    float y_row_1 = r13*(P.x()-N.x()) + r23*(P.y() - N.y()) + r33*(P.z() - N.z());

    float x = H.x() - -focal_length * (x_row_0 / x_row_1);
    float y = H.y() - -focal_length * (y_row_0 / y_row_1);
    QVector2D projected_point_without_z = QVector2D(x, y);

    return move_into_plane(projected_point_without_z);
}

std::vector<QVector3D> Camera::project(Cube cube, int algorithm) {
    std::vector<QVector3D> arr;
    for (auto vertex : cube.get_cube_points()) {
        if (algorithm == 0) {
            arr.push_back(affine_projection(vertex));
        } else if (algorithm == 1) {
            arr.push_back(homogeneous_projection(vertex));
        } else {
            arr.push_back(explicit_projection(vertex));
        }
    }
    return arr;
}
