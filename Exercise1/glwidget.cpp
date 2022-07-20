
// (c) Nico Brügel, 2021

#include "glwidget.h"
#include <QtGui>

#if defined(__APPLE__)
// we're on macOS and according to their documentation Apple hates developers
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
// windows, even if it's case insensitive QT-Create might generate a warning
#include <gl/GL.h>
#include <gl/GLU.h>
#else
// hopefully on linux
// If can't be found, ensure that the following is installed:
// libglu1-mesa-dev and/or mesa-common-dev
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <QMouseEvent>
#include <QFileDialog>
#include <QMessageBox>

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <utility>

#include "mainwindow.h"
#include "GLConvenience.h"
#include "QtConvenience.h"

#include <Eigen/Eigenvalues>

using namespace std;
#define PI 3.14159265

GLWidget::GLWidget(QWidget* parent)
    : QOpenGLWidget(parent),
    _pointSize(1)
{
    setMouseTracking(true);

    // axes cross
    _axesLines.push_back(make_pair(QVector3D(0.0, 0.0, 0.0), QColor(1.0, 0.0, 0.0)));
    _axesLines.push_back(make_pair(QVector3D(1.0, 0.0, 0.0), QColor(1.0, 0.0, 0.0)));
    _axesLines.push_back(make_pair(QVector3D(0.0, 0.0, 0.0), QColor(0.0, 1.0, 0.0)));
    _axesLines.push_back(make_pair(QVector3D(0.0, 1.0, 0.0), QColor(0.0, 1.0, 0.0)));
    _axesLines.push_back(make_pair(QVector3D(0.0, 0.0, 0.0), QColor(0.0, 0.0, 1.0)));
    _axesLines.push_back(make_pair(QVector3D(0.0, 0.0, 1.0), QColor(0.0, 0.0, 1.0)));
}

GLWidget::~GLWidget()
{
    this->cleanup();
}

void GLWidget::cleanup()
{
  makeCurrent();
 // _vertexBuffer.destroy();
  _shaders.reset();
  doneCurrent();
}

void GLWidget::initializeGL()
{
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);

  initializeOpenGLFunctions();
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(0.4f,0.4f,0.4f,1);     // screen background color

  // the world is still for now
  _worldMatrix.setToIdentity();

  // create shaders and map attributes
  initShaders();

  // create array container and load points into buffer
  createContainers();
}

void GLWidget::paintGL()
{
    // ensure GL flags
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE); //required for gl_PointSize

    // create shaders and map attributes
    initShaders();

    // create array container and load points into buffer
   // createContainers();
    // create array container and load points into buffer
    const QVector<float>& pointsData =pointcloud.getData();
    //const QVector<float>& pointsData = source_pointcloud.getData();
    //const QVector<float>& pointsData2 = target_pointcloud.getData();
    if(!_vao.isCreated()) _vao.create();
    QOpenGLVertexArrayObject::Binder vaoBinder(&_vao);
    if(!_vertexBuffer.isCreated()) _vertexBuffer.create();
    _vertexBuffer.bind();
    _vertexBuffer.allocate(pointsData.constData(), pointsData.size() * sizeof(GLfloat));
    //_vertexBuffer.allocate(pointsData3.constData(), pointsData3.size() * sizeof(GLfloat));
    QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
    f->glEnableVertexAttribArray(0);
    f->glEnableVertexAttribArray(1);
    f->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat) + sizeof(GLfloat), nullptr);
    f->glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat) + sizeof(GLfloat), reinterpret_cast<void *>(3*sizeof(GLfloat)));
    _vertexBuffer.release();

    // set camera
    setupRenderingCamera();

    // draw points cloud
    drawPointCloud();
    draw_source_pointcloud();
    draw_target_pointcloud();
    drawFrameAxis();

    // Assignement 1, Part 1
    // Draw here your objects as in drawFrameAxis();
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //glBegin(GL_LINES);
    Cube cube1 = Cube(QVector3D(5, 0, 4), QVector3D(45, 45, 0), 0.5);
    Cube cube2 = Cube(QVector3D(3, 0, 4), QVector3D(0, 0, 0), 0.5);
    Cube cube3 = Cube(QVector3D(4, 0, 20), QVector3D(0, 0, 0), 0.5);

    if(_cubes_on){
        draw_cube(cube1, QColor(0.0, 1.0, 0.0));
        draw_cube(cube2, QColor(0.0, 1.0, 0.0));
        draw_cube(cube3, QColor(0.0, 1.0, 0.0));
    }

    // Assignement 1, Part 2
    // Draw here your perspective camera model
    // Assignement 1, Part 3
    // Draw here the perspective projection
    Camera camera1 = Camera(QVector3D(_c1_x, _c1_y, _c1_z), QVector3D(_c1_rotation_x, _c1_rotation_y, _c1_rotation_z), _focal_length);
    Camera camera2 = Camera(QVector3D(0, 0, 0), QVector3D(0, 0, 0), 2);

    std::vector<QVector3D> c1_projected_1;
    std::vector<QVector3D> c1_projected_2;
    std::vector<QVector3D> c1_projected_3;

    std::vector<QVector3D> c2_projected_1;
    std::vector<QVector3D> c2_projected_2;
    std::vector<QVector3D> c2_projected_3;

    if (_rays_on) {
        if (_camera1_on) {
            draw_rays(camera1, cube1);
            draw_rays(camera1, cube2);
            draw_rays(camera1, cube3);
        }

        if (_camera2_on) {
            draw_rays(camera2, cube1);
            draw_rays(camera2, cube2);
            draw_rays(camera2, cube3);
        }
    }

    if (_camera1_on) {
        if (_projection_on) {
            std::vector<QVector3D> c1_projected_1 = camera1.project(cube1, _projection_algorithm);
            std::vector<QVector3D> c1_projected_2 = camera1.project(cube2, _projection_algorithm);
            std::vector<QVector3D> c1_projected_3 = camera1.project(cube3, _projection_algorithm);

            draw_projection(c1_projected_1);
            draw_projection(c1_projected_2);
            draw_projection(c1_projected_3);
        }
        draw_camera_model(camera1, 1, 2.5);
    }

    if (_camera2_on) {
        if (_projection_on) {
            std::vector<QVector3D> c2_projected_1 = camera2.project(cube1, _projection_algorithm);
            std::vector<QVector3D> c2_projected_2 = camera2.project(cube2, _projection_algorithm);
            std::vector<QVector3D> c2_projected_3 = camera2.project(cube3, _projection_algorithm);

            draw_projection(c2_projected_1);
            draw_projection(c2_projected_2);
            draw_projection(c2_projected_3);
        }
        draw_camera_model(camera2, 1, 5);
    }

    // Assignment 2, Part 2
    if (_camera1_on && _camera2_on) {
        perspective_reconstruction(camera2, camera1, cube1);
        perspective_reconstruction(camera2, camera1, cube2);
        perspective_reconstruction(camera2, camera1, cube3);
    }

    // Assignment 3
    // 1.1) construction of kd-tree
    // is done if data is laded
    if (draw_kd_tree_partitions) {
        std::vector<KdNode> nodes = kdtree.get_all_level_nodes(2);
        std::vector<QVector3D> points;

        for (auto node : nodes) {
            points.push_back(node.point);
        }
        draw_points(points, 15, QColor(1.0, 1.0, 0.0));

        for (auto node : nodes) {
            draw_plane(node.point, node.direction, 5);
        }
    }

    if(draw_oct_tree_partitions) {
        std::vector<OctNode> list_of_octnodes =  octtree.get_all_level_nodes(1);

        for (auto node : list_of_octnodes) {
            draw_cube_2(Cube(node.bottomRightBack ,node.topLeftFront), QColor(1.0, 1.0, 1.0));
        }
        //draw_points(octtree.test);
    }

    // Assignment 4 PCA
    if (source_pointcloud.getCount() != 0) {
        draw_pca_axes(&source_pointcloud);
    }

    if (target_pointcloud.getCount() != 0) {
        //target_pointcloud.test = true;
        draw_pca_axes(&target_pointcloud);
    }

    update_source_pointcloud();
}

void GLWidget::update_source_pointcloud() {
    float x = _c1_rotation_x * M_PI / 180;
    float y = _c1_rotation_y * M_PI / 180;
    float z = _c1_rotation_z * M_PI / 180;

    QMatrix4x4 rotation_matrix_x = QMatrix4x4(1, 0, 0, 0, 0, cos(x), -sin(x), 0, 0, sin(x), cos(x), 0, 0,0,0,1);
    QMatrix4x4 rotation_matrix_y = QMatrix4x4(cos(y), 0, sin(y), 0, 0, 1, 0, 0, -sin(y), 0, cos(y), 0, 0,0,0,1 );
    QMatrix4x4 rotation_matrix_z = QMatrix4x4(cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0,0,0,0,1);
    QMatrix4x4 rotation_matrix = rotation_matrix_z * rotation_matrix_y * rotation_matrix_x;
    QVector3D translation = QVector3D(_c1_x, _c1_y, _c1_z);
    //rotation_and_translation_matrix.setColumn(3, QVector4D(-0.3, 0, 0, 1));

    std::vector<QVector3D> new_points;
    for (auto point: target_pointcloud.getPoints()) {
        //QVector3D new_point = rotation_and_translation_matrix * point;
        QVector3D new_point = translation + target_pointcloud.centroid + rotation_matrix * (point - target_pointcloud.centroid);
        new_points.push_back(new_point);
    }
    source_pointcloud.set_points(new_points);
}

void GLWidget::draw_pca_axes(PointCloud* pc) {
    // Get points and calculate centroid Point
    float cx = 0, cy = 0, cz = 0;
    for (auto point: pc->getPoints()) {
        cx += point.x();
        cy += point.y();
        cz += point.z();
    }
    int m = pc->getPoints().size();
    QVector3D centroid = QVector3D(cx/m, cy/m, cz/m);
    pc->centroid = centroid;
    draw_point(centroid);

    // Set up the big matrix
    std::vector<float> px, py, pz;
    for (auto point: pc->getPoints()) {
        px.push_back(point.x() - centroid.x());
        py.push_back(point.y() - centroid.y());
        pz.push_back(point.z() - centroid.z());
    }

    // Calculate M matrix, which is symmetric!
    QMatrix4x4 M = QMatrix4x4(multiply_rows(px, px), multiply_rows(px, py), multiply_rows(px, pz), 0,
                              multiply_rows(py, px), multiply_rows(py, py), multiply_rows(py, pz), 0,
                              multiply_rows(pz, px), multiply_rows(pz, py), multiply_rows(pz, pz), 0,
                              0, 0, 0, 1);

    // Option A: get eigenvalues and eigenvectors
    std::vector<float> vals = eigvals(M);
    std::vector<QVector3D> vecs = eigvecs(M, vals);

    QMatrix4x4 matrix = QMatrix4x4();
    matrix.setColumn(0, vecs[0].toVector4D());
    matrix.setColumn(1, vecs[1].toVector4D());
    matrix.setColumn(2, vecs[2].toVector4D());
    matrix.setColumn(3, QVector4D(0,0,0,1));
    // Ende Option A

    // Option B:
//    Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
//    Eigen::Matrix3f covmat;
//    covmat << multiply_rows(px, px), multiply_rows(px, py), multiply_rows(px, pz),
//         multiply_rows(py, px), multiply_rows(py, py), multiply_rows(py, pz),
//         multiply_rows(pz, px), multiply_rows(pz, py), multiply_rows(pz, pz);
//    eigensolver.compute(covmat);
//    Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
//    Eigen::MatrixXf eigen_vectors = eigensolver.eigenvectors().real();
//    std::vector<std::tuple<float, Eigen::VectorXf>> eigen_vectors_and_values;

//    for(int i=0; i<eigen_values.size(); i++){
//        std::tuple<float, Eigen::VectorXf> vec_and_val(eigen_values[i], eigen_vectors.row(i));
//        eigen_vectors_and_values.push_back(vec_and_val);
//    }
//    std::sort(eigen_vectors_and_values.begin(), eigen_vectors_and_values.end(),
//        [&](const std::tuple<float, Eigen::VectorXf>& a, const std::tuple<float, Eigen::VectorXf>& b) -> bool{
//            return std::get<0>(a) <= std::get<0>(b);
//    });
//    int index = 0;
//    for(auto const vect : eigen_vectors_and_values){
//        eigen_values(index) = std::get<0>(vect);
//        eigen_vectors.row(index) = std::get<1>(vect);
//        index++;
//    }

//    QMatrix4x4 matrix = QMatrix4x4();
//    matrix.setColumn(0, QVector4D(eigen_vectors.row(0).x(), eigen_vectors.row(0).y(), eigen_vectors.row(0).z(), 0));
//    matrix.setColumn(1, QVector4D(eigen_vectors.row(1).x(), eigen_vectors.row(1).y(), eigen_vectors.row(1).z(), 0));
//    matrix.setColumn(2, QVector4D(eigen_vectors.row(2).x(), eigen_vectors.row(2).y(), eigen_vectors.row(2).z(), 0));
//    matrix.setColumn(3, QVector4D(0,0,0,1));
//    std::vector<QVector3D> vecs = {QVector3D(eigen_vectors.row(0).x(), eigen_vectors.row(0).y(), eigen_vectors.row(0).z()),
//                                  QVector3D(eigen_vectors.row(1).x(), eigen_vectors.row(1).y(), eigen_vectors.row(1).z()),
//                                  QVector3D(eigen_vectors.row(2).x(), eigen_vectors.row(2).y(), eigen_vectors.row(2).z())};
    // Ende Option B

    pc->eigvecs_matrix = matrix;
    // draw the axes
    draw_axes(vecs, centroid, 0.1);
}

std::vector<float> GLWidget::eigvals(QMatrix4x4 matrix) {
    std::vector<float> vals;
    float a = matrix.row(0).x();
    float b = matrix.row(0).y();
    float c = matrix.row(0).z();
    float d = matrix.row(1).x();
    float e = matrix.row(1).y();
    float f = matrix.row(1).z();
    float g = matrix.row(2).x();
    float h = matrix.row(2).y();
    float i = matrix.row(2).z();

    float trace = a + e + i;
    float determinant = (a*e*i - a*f*h) + (b*f*g - b*d*i) + (c*d*h - c*e*g);

    std::vector<float> px = {matrix.row(0).x(), matrix.row(0).y(), matrix.row(0).z()};
    std::vector<float> py = {matrix.row(1).x(), matrix.row(1).y(), matrix.row(1).z()};
    std::vector<float> pz = {matrix.row(2).x(), matrix.row(2).y(), matrix.row(2).z()};

    QMatrix4x4 second_matrix = QMatrix4x4(multiply_rows(px, px), multiply_rows(px, py), multiply_rows(px, pz), 0,
                              multiply_rows(py, px), multiply_rows(py, py), multiply_rows(py, pz), 0,
                              multiply_rows(pz, px), multiply_rows(pz, py), multiply_rows(pz, pz), 0,
                              0, 0, 0, 1);
    float second_trace = second_matrix.row(0).x() + second_matrix.row(1).y() + second_matrix.row(2).z();
    std::vector<float> test = solve_cubic_equation(-1, trace, -0.5*(trace*trace - second_trace), determinant);
    return solve_cubic_equation(-1, trace, -0.5*(trace*trace - second_trace), determinant);
}

std::vector<float> GLWidget::solve_cubic_equation(float a, float b, float c, float d) {
    std::vector<float> roots;

    float p = (b*b - 3*a*c);
    float q = (9.*a*b*c - 2*(b*b*b) - 27*(a*a)*d);
    //float n = (27*(p*p*p) / (q*q));
    float theta = acos(q / ((2*p)*sqrt(p)));
    float root1 = (-b + 2*cos(theta / 3)*sqrt(p)) / (3*a);
    float root2 = (-b + 2*cos((theta / 3) + (120. * (PI / 180)))*sqrt(p)) / (3*a);
    float root3 = (-b + 2*cos((theta / 3) + (240. * (PI / 180)))*sqrt(p)) / (3*a);

    roots.push_back(root1);
    roots.push_back(root2);
    roots.push_back(root3);
    std::sort(roots.begin(), roots.end(), std::greater{});
    return roots;
}

std::vector<QVector3D> GLWidget::eigvecs(QMatrix4x4 matrix, std::vector<float> vals) {
    std::vector<QVector3D> vecs;

    float a = matrix.row(0).x();
    float b = matrix.row(0).y();
    float c = matrix.row(0).z();
    float d = matrix.row(1).x();
    float e = matrix.row(1).y();
    float f = matrix.row(1).z();
    float g = matrix.row(2).x();
    float h = matrix.row(2).y();
    float i = matrix.row(2).z();

    for (auto val: vals) {
        std::vector<float> row0 = {a - val, b, c};
        std::vector<float> row1 = {d, e - val, f};
        std::vector<float> row2 = {g, h, i - val};

        // Kreuzprodukt
        QVector3D r0xr1 = QVector3D(row0[1] * row1[2] - row0[2] * row1[1],
                                    row0[2] * row1[0] - row0[0] * row1[2],
                                    row0[0] * row1[1] - row0[1] * row1[0]);
        QVector3D r0xr2 = QVector3D(row0[1] * row2[2] - row0[2] * row2[1],
                                    row0[2] * row2[0] - row0[0] * row2[2],
                                    row0[0] * row2[1] - row0[1] * row2[0]);
        QVector3D r1xr2 = QVector3D(row1[1] * row2[2] - row1[2] * row2[1],
                                    row1[2] * row2[0] - row1[0] * row2[2],
                                    row1[0] * row2[1] - row1[1] * row2[0]);
        float d0 = r0xr1[0]*r0xr1[0] + r0xr1[1]*r0xr1[1] + r0xr1[2]*r0xr1[2];
        float d1 = r0xr2[0]*r0xr2[0] + r0xr2[1]*r0xr2[1] + r0xr2[2]*r0xr2[2];
        float d2 = r1xr2[0]*r1xr2[0] + r1xr2[1]*r1xr2[1] + r1xr2[2]*r1xr2[2];

        float index_max = 0;
        float d_max = d0;

        if (d1 > d_max) {
            d_max = d1;
            index_max = 1;
        }
        if (d2 > d_max) {
            d_max = d2;
            index_max = 2;
        }

        if (index_max == 0) {
            //vecs.push_back(r0xr1 / r0xr1[2]);
            vecs.push_back(r0xr1 / sqrt(d0));
        }
        if (index_max == 1) {
            //vecs.push_back(r0xr2 / r0xr2[2]);
            vecs.push_back(r0xr2 / sqrt(d1));
        }
        if (index_max == 2) {
            //vecs.push_back(r1xr2 / r1xr2[2]);
            vecs.push_back(r1xr2 / sqrt(d2));
        }
    }
    return vecs;
}

float GLWidget::multiply_rows(std::vector<float> row1, std::vector<float> row2) {
    float sum = 0;
    for (int i = 0; i < row1.size(); i++) {
        sum += row1[i] * row2[i];
    }
    return sum;
}

void GLWidget::draw_point(QVector3D point) {
    glPointSize(10);
    glBegin(GL_POINTS);
    const auto translated = _projectionMatrix * _cameraMatrix * _worldMatrix ^ point;
    glColor3f(QColor(1.0, 1.0, 0.0));
    glVertex3f(translated);
    glEnd();
}

void GLWidget::draw_axes(std::vector<QVector3D> axes, QVector3D centroid, float scale) {
    //axes
    QVector3D x_axis = axes[0];
    QVector3D y_axis = axes[1];
    QVector3D z_axis = axes[2];

    //x_axis = x_axis + (centroid - x_axis) * scale;
    //y_axis = y_axis + (centroid - y_axis) * scale;
    //z_axis = z_axis + (centroid - z_axis) * scale;
    x_axis = centroid + x_axis * scale;
    y_axis = centroid + y_axis * scale;
    z_axis = centroid + z_axis * scale;

    std::vector<std::pair<QVector3D, QColor>> axes_lines;
    axes_lines.push_back(make_pair(centroid, QColor(1.0, 0.0, 0.0)));
    axes_lines.push_back(make_pair(x_axis, QColor(1.0, 0.0, 0.0)));

    axes_lines.push_back(make_pair(centroid, QColor(0.0, 1.0, 0.0)));
    axes_lines.push_back(make_pair(y_axis, QColor(0.0, 1.0, 0.0)));

    axes_lines.push_back(make_pair(centroid, QColor(0.0, 0.0, 1.0)));
    axes_lines.push_back(make_pair(z_axis, QColor(0.0, 0.0, 1.0)));

    // Draw axes
    glBegin(GL_LINES);
    for (auto vertex : axes_lines) {
      const auto translated = _projectionMatrix * _cameraMatrix * _worldMatrix ^ vertex.first;
      glColor3f(vertex.second);
      glVertex3f(translated);
    }
    glEnd();
}

void GLWidget::draw_plane(QVector3D point, int direction, double size) {
    QVector4D direction1, direction2, direction3, direction4;
    float x = point.x(), y = point.y(), z = point.z();
    float scaled_size = size / 100;

    glLineWidth(1);
    glBegin(GL_QUADS);

    if (direction == 0) { // 0 = X
        glColor4f(QColor(1.0, 0.0, 0.0), 0.4);
        direction1 = QVector4D(x, y + scaled_size, z + scaled_size, 0);
        direction2 = QVector4D(x, y + scaled_size, z - scaled_size, 0);
        direction3 = QVector4D(x, y - scaled_size, z - scaled_size, 0);
        direction4 = QVector4D(x, y - scaled_size, z + scaled_size, 0);

    } else if (direction == 1) { // 1 = Y
        glColor4f(QColor(0.0, 1.0, 0.0), 0.4);
        direction1 = QVector4D(x + scaled_size, y, z + scaled_size, 0);
        direction2 = QVector4D(x + scaled_size, y, z - scaled_size, 0);
        direction3 = QVector4D(x - scaled_size, y, z - scaled_size, 0);
        direction4 = QVector4D(x - scaled_size, y, z + scaled_size, 0);
    } else if (direction == 2) { // 2 = Z
        glColor4f(QColor(0.0, 0.0, 1.0), 0.4);
        direction1 = QVector4D(x + scaled_size, y + scaled_size, z, 1);
        direction2 = QVector4D(x + scaled_size, y - scaled_size, z, 1);
        direction3 = QVector4D(x - scaled_size, y - scaled_size, z, 1);
        direction4 = QVector4D(x - scaled_size, y + scaled_size, z, 1);
    }

    QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
    const auto direction1_p = _projectionMatrix * mvMatrix ^ direction1;
    const auto direction2_p = _projectionMatrix * mvMatrix ^ direction2;
    const auto direction3_p = _projectionMatrix * mvMatrix ^ direction3;
    const auto direction4_p = _projectionMatrix * mvMatrix ^ direction4;

    glVertex3f(direction1_p);
    glVertex3f(direction2_p);
    glVertex3f(direction3_p);
    glVertex3f(direction4_p);

    glEnd();
}

void GLWidget::draw_points(std::vector<QVector3D> points, int size, QColor color) {
    glPointSize(size);
    glBegin(GL_POINTS);
    for (auto point : points) {
        const auto translated = _projectionMatrix * _cameraMatrix * _worldMatrix ^ point;
        glColor3f(QColor(color));
        glVertex3f(translated);
    }
    glEnd();
}

void GLWidget::perspective_reconstruction(Camera left_camera, Camera right_camera, Cube cube) {
    float b = right_camera.getPosition().x();
    float f = left_camera.getF();

    std::vector<QVector3D> left_projection = left_camera.project(cube, _projection_algorithm);
    std::vector<QVector3D> right_projection = right_camera.project(cube, _projection_algorithm);

    std::vector<QVector3D> points;

    QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
    mvMatrix.scale(0.05f);

    glPointSize(10);
    glBegin(GL_POINTS);
    for (int i = 0; i < left_projection.size(); i++) {
        QVector3D left_point = left_projection.at(i) - left_camera.getImagePrinciplePoint();
        QVector3D right_point = right_projection.at(i) - right_camera.getImagePrinciplePoint();

        float x_parallax = left_point.x() - right_point.x() ;

        float z = f * b / x_parallax; // before -f
        float y = -z * left_point.y() / -f; // before f
        float x = -z * left_point.x() / -f; // before f

        QVector3D reconstructed_point = _projectionMatrix * mvMatrix ^ QVector3D(x, y, z);
        points.push_back(reconstructed_point);

        glColor3f(QColor(1.0, 0, 0));
        glVertex3f(reconstructed_point);
    }
    glEnd();

    glBegin(GL_LINES);
    for (int i = 0; i < points.size(); i++) {
        QVector3D reconstructed_point = points.at(i);
        glColor3f(QColor(1.0, 0, 0));
        glVertex3f(reconstructed_point);
    }
    glEnd();
}

void GLWidget::draw_projection(std::vector<QVector3D> projected_points) {
    QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
    mvMatrix.scale(0.05f);

    glBegin(GL_LINES);
    for (auto point : projected_points) {
        QVector3D to_draw = _projectionMatrix * mvMatrix ^ point;
        glColor3f(QColor(0.0, 1.0, 0.0));
        glVertex3f(to_draw);
    }
    glEnd();

    glPointSize(10);
    glBegin(GL_POINTS);
    for (auto point : projected_points) {
        QVector3D to_draw = _projectionMatrix * mvMatrix ^ point;
        glColor3f(QColor(0.0, 1.0, 0.0));
        glVertex3f(to_draw);
    }
    glEnd();
}

void GLWidget::draw_rays(Camera camera, Cube cube) {
    QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
    mvMatrix.scale(0.05f);

    glBegin(GL_LINES);
    for (auto point : cube.get_cube_points()) {
      auto cube_point = _projectionMatrix * mvMatrix ^ point;
      auto camera_point = _projectionMatrix * mvMatrix ^ camera.getPosition();
      glColor3f(QColor(1.0, 1.0, 1.0));
      glVertex3f(cube_point);
      glVertex3f(camera_point);
    }
    glEnd();
}

void GLWidget::draw_camera_model(Camera camera, float height, float width){
    float hX = camera.getImagePrinciplePoint().x();
    float hY = camera.getImagePrinciplePoint().y();
    float hZ = camera.getImagePrinciplePoint().z();

    QVector3D camera_position = camera.getPosition();
    QVector4D projection_center = QVector4D(camera_position.x(), camera_position.y(), camera_position.z(), 0);
    QVector4D image_principle_point = QVector4D(hX, hY, hZ, 0);
    QMatrix4x4 transformation_matrix = camera.getTransformationMatrix();
    QMatrix4x4 rotation_matrix = camera.getRotationMatrix();

    // TODO move logic to camera class
    //axes
    QVector4D x_axes = QVector4D(1, 0, 0, 1);
    x_axes = transformation_matrix * x_axes;

    QVector4D y_axes = QVector4D(0, 1, 0, 1);
    y_axes = transformation_matrix * y_axes;

    QVector4D z_axes = QVector4D(0, 0, 1, 1);
    z_axes = transformation_matrix * z_axes;

    std::vector<std::pair<QVector3D, QColor>> axes_lines;
    axes_lines.push_back(make_pair(camera_position, QColor(1.0, 0.0, 0.0)));
    axes_lines.push_back(make_pair(x_axes.toVector3D(), QColor(1.0, 0.0, 0.0)));

    axes_lines.push_back(make_pair(camera_position, QColor(0.0, 1.0, 0.0)));
    axes_lines.push_back(make_pair(y_axes.toVector3D(), QColor(0.0, 1.0, 0.0)));

    axes_lines.push_back(make_pair(camera_position, QColor(0.0, 0.0, 1.0)));
    axes_lines.push_back(make_pair(z_axes.toVector3D(), QColor(0.0, 0.0, 1.0)));

    QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
    mvMatrix.scale(0.05f); // make it small

    // Draw axes
    glBegin(GL_LINES);
    for (auto vertex : axes_lines) {
      const auto translated = _projectionMatrix * mvMatrix ^ vertex.first;
      glColor3f(vertex.second);
      glVertex3f(translated);
    }
    glEnd();

    // Draw projection center point and image principle point
    glPointSize(15);
    glBegin(GL_POINTS);
        auto projection_center_translated = _projectionMatrix * mvMatrix ^ camera.getPosition();
        glColor3f(QColor(1.0, 1.0, 0.0));
        glVertex3f(projection_center_translated);
        auto image_principle_point_translated = _projectionMatrix * mvMatrix ^ camera.getImagePrinciplePoint();
        glVertex3f(image_principle_point_translated);
    glEnd();
    glPointSize(1);

    // Draw line from image principle point to projection center point
    glLineWidth(5);
    glBegin(GL_LINES);
        glVertex3f(image_principle_point_translated);
        QVector3D image_principle_vector_normalized =  (projection_center.toVector3D()-image_principle_point.toVector3D()).normalized();
        auto image_principle_vector_translated = _projectionMatrix * mvMatrix ^ (image_principle_point.toVector3D() + (image_principle_vector_normalized * 0.33));
        glVertex3f(image_principle_vector_translated);
   glEnd();

   // Draw image plane relative to camera pose
   glLineWidth(1);
   glBegin(GL_QUADS);
        glColor4f(QColor(1.0, 1.0, 0.0),0.5);

        // TODO move logic to camera class
        QVector4D direction1 = QVector4D(width / 2, height / 2, 0, 1);
        QVector4D direction2 = QVector4D(width / 2, -height / 2, 0, 1);
        QVector4D direction3 = QVector4D(-width / 2, -height / 2, 0, 1);
        QVector4D direction4 = QVector4D(-width / 2, height / 2, 0, 1);

        auto imagespace_point_translated = _projectionMatrix * mvMatrix ^ (image_principle_point + rotation_matrix * direction1).toVector3D();
        glVertex3f(imagespace_point_translated);
        imagespace_point_translated = _projectionMatrix * mvMatrix ^ (image_principle_point + rotation_matrix * direction2).toVector3D();
        glVertex3f(imagespace_point_translated);
        imagespace_point_translated = _projectionMatrix * mvMatrix ^ (image_principle_point + rotation_matrix * direction3).toVector3D();
        glVertex3f(imagespace_point_translated);
        imagespace_point_translated = _projectionMatrix * mvMatrix ^ (image_principle_point + rotation_matrix * direction4).toVector3D();
        glVertex3f(imagespace_point_translated);
   glEnd();
}

void GLWidget::draw_cube(Cube cube, QColor color)
{
  glBegin(GL_LINES);
  QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
  mvMatrix.scale(0.05f); // make it small
  for (auto vertex : cube.get_cube_points()) {
    const auto translated = _projectionMatrix * mvMatrix ^ vertex;
    glColor3f(color);
    glVertex3f(translated);
  }
  glEnd();
}

void GLWidget::draw_cube_2(Cube cube, QColor color)
{
  glBegin(GL_LINES);
  QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
  mvMatrix.scale(0.05f); // make it small
  for (auto vertex : cube.get_cube_from_2_points()) {
    const auto translated = _projectionMatrix * mvMatrix ^ vertex;
    glColor3f(color);
    glVertex3f(translated);
  }
  glEnd();

  bool debugpoints = false;
  if(debugpoints){
      glPointSize(15);
      glBegin(GL_POINTS);
          auto translated = _projectionMatrix * _cameraMatrix * _worldMatrix ^ cube.position;
          glColor3f(QColor(1.0, 1.0, 0.0));
          glVertex3f(translated);
      glEnd();
      glPointSize(15);
      glBegin(GL_POINTS);
          translated = _projectionMatrix * _cameraMatrix * _worldMatrix ^ cube.position_2;
          glColor3f(QColor(1.0, 0.0, 0.0));
          glVertex3f(translated);
      glEnd();
  }

}

void GLWidget::drawFrameAxis()
{
  glBegin(GL_LINES);
  QMatrix4x4 mvMatrix = _cameraMatrix * _worldMatrix;
  mvMatrix.scale(0.05f); // make it small
  for (auto vertex : _axesLines) {
    const auto translated = _projectionMatrix * mvMatrix ^ vertex.first;
    glColor3f(vertex.second);
    glVertex3f(translated);
  }
  glEnd();
}

void GLWidget::resizeGL(int w, int h)
{
  _projectionMatrix.setToIdentity();
  _projectionMatrix.perspective(70.0f, GLfloat(w) / GLfloat(h), 0.01f, 100.0f);
}

void GLWidget::wheelEvent(QWheelEvent* event)
{
  if (event->angleDelta().y() > 0) {
    _renderingCamera->forward();
  } else {
    _renderingCamera->backward();
  }
}

void GLWidget::keyPressEvent(QKeyEvent * event)
{
    switch ( event->key() )
    {
      case Qt::Key_Escape:
        QApplication::instance()->quit();
        break;

      case Qt::Key_Left:
      case Qt::Key_A:
        _renderingCamera->left();
        break;

      case Qt::Key_Right:
      case Qt::Key_D:
        _renderingCamera->right();
        break;

      case Qt::Key_Up:
      case Qt::Key_W:
        _renderingCamera->forward();
        break;

      case Qt::Key_Down:
      case Qt::Key_S:
        _renderingCamera->backward();
        break;

      case Qt::Key_Space:
      case Qt::Key_Q:
        _renderingCamera->up();
        break;

      case Qt::Key_C:
      case Qt::Key_Z:
        _renderingCamera->down();
        break;

      default:
        QWidget::keyPressEvent(event);
    }
    update();
}


void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    const QPoint pos = event->pos() - _prevMousePosition;

    _prevMousePosition = event->pos();

    if (event->buttons() & Qt::LeftButton)
    {
        _renderingCamera->rotate(X_Pressed?0:pos.y(), Y_Pressed?0:pos.x(), 0);
    }
    else if ( event->buttons() & Qt::RightButton)
    {
        if (pos.x() < 0) _renderingCamera->right();
        if (pos.x() > 0) _renderingCamera->left();
        if (pos.y() < 0) _renderingCamera->down();
        if (pos.y() > 0) _renderingCamera->up();
    }
}

void GLWidget::attachCamera(QSharedPointer<RenderingCamera> camera)
{
  if (_renderingCamera)
  {
    disconnect(_renderingCamera.data(), &RenderingCamera::changed, this, &GLWidget::onCameraChanged);
  }
  _renderingCamera = camera;
  connect(camera.data(), &RenderingCamera::changed, this, &GLWidget::onCameraChanged);
}


void GLWidget::onCameraChanged(const RenderingCameraState&)
{
  update();
}

void GLWidget::setPointSize(size_t size)
{
  assert(size > 0);
  _pointSize = static_cast<float>(size);
  update();
}

void GLWidget::openFileDialog()
{
    const QString filePath = QFileDialog::getOpenFileName(this, tr("Open PLY file"), "../data", tr("PLY Files (*.ply)"));
     if (!filePath.isEmpty())
     {
         std::cout << filePath.toStdString() << std::endl;
         pointcloud.loadPLY(filePath);

         // kd Tree
         kdtree = KdTree(pointcloud._X, pointcloud._Y, pointcloud._Z);
         kdtree.construct_balanced_tree();
         //draw_kd_tree_partitions = true;

         //Quad Tree
         octtree = OctTree(pointcloud._X, pointcloud._Y, pointcloud._Z);
         octtree.construct_balanced_tree();
         octtree_root = octtree.root;

         update();
     }
}

// Open target point cloud must be opened first!
void GLWidget::open_source_pointcloud() {
    const QString filePath = QFileDialog::getOpenFileName(this, tr("Open source pointcloud"), "../data", tr("PLY Files (*.ply)"));
    if (!filePath.isEmpty()) {
        std::cout << filePath.toStdString() << std::endl;
        source_pointcloud.loadPLY(filePath);

        float x = _c1_rotation_x * M_PI / 180;
        float y = _c1_rotation_y * M_PI / 180;
        float z = _c1_rotation_z * M_PI / 180;

        QMatrix4x4 rotation_matrix_x = QMatrix4x4(1, 0, 0, 0, 0, cos(x), -sin(x), 0, 0, sin(x), cos(x), 0, 0,0,0,1);
        QMatrix4x4 rotation_matrix_y = QMatrix4x4(cos(y), 0, sin(y), 0, 0, 1, 0, 0, -sin(y), 0, cos(y), 0, 0,0,0,1 );
        QMatrix4x4 rotation_matrix_z = QMatrix4x4(cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0,0,0,0,1);
        QMatrix4x4 rotation_matrix = rotation_matrix_z * rotation_matrix_y * rotation_matrix_x;
        QVector3D translation = QVector3D(_c1_x, _c1_y, _c1_z);
        //rotation_and_translation_matrix.setColumn(3, QVector4D(-0.3, 0, 0, 1));

        std::vector<QVector3D> new_points;
        for (auto point: target_pointcloud.getPoints()) {
            //QVector3D new_point = rotation_and_translation_matrix * point;
            QVector3D new_point = translation + target_pointcloud.centroid + rotation_matrix * (point - target_pointcloud.centroid);
            new_points.push_back(new_point);
        }
        source_pointcloud.set_points(new_points);
        update();
    }
}

void GLWidget::open_target_pointcloud() {
    const QString filePath = QFileDialog::getOpenFileName(this, tr("Open target pointcloud"), "../data", tr("PLY Files (*.ply)"));
    if (!filePath.isEmpty()) {
        std::cout << filePath.toStdString() << std::endl;
        target_pointcloud.loadPLY(filePath);
        update();
    }
}

void GLWidget::align_point_clouds() {
    // R * A = B
    // R = B * A^-1 = B * A^T
    // A = source matrix
    // B = target matrix

    // Get rotation matrix with smalles angle
    float trace = 0;
    QMatrix4x4 R;

    QMatrix4x4 source_matrix = source_pointcloud.eigvecs_matrix;
    QMatrix4x4 R1 = target_pointcloud.eigvecs_matrix * source_matrix.transposed();
    float trace1 = calculate_angle_from_trace(R1);
    trace = trace1;
    R = R1;

//    // Invert second and third eigenvector
//    QMatrix4x4 source_matrix_2;
//    source_matrix_2.setColumn(0, source_matrix.column(0));
//    source_matrix_2.setColumn(1, source_matrix.column(1) * -1);
//    source_matrix_2.setColumn(2, source_matrix.column(2) * -1);
//    QMatrix4x4 R2 = target_pointcloud.eigvecs_matrix * source_matrix_2.transposed();
//    float trace2 = calculate_angle_from_trace(R2);
//    if (trace2 < trace) {
//        trace = trace2;
//        R = R2;
//    }

//    // Invert first and third eigenvector
//    QMatrix4x4 source_matrix_3;
//    source_matrix_3.setColumn(0, source_matrix.column(0) * -1);
//    source_matrix_3.setColumn(1, source_matrix.column(1));
//    source_matrix_3.setColumn(2, source_matrix.column(2) * -1);
//    QMatrix4x4 R3 = target_pointcloud.eigvecs_matrix * source_matrix_3.transposed();
//    float trace3 = calculate_angle_from_trace(R3);
//    if (trace3 < trace) {
//        trace = trace3;
//        R = R3;
//    }

//    // Invert first and second eigenvector
//    QMatrix4x4 source_matrix_4;
//    source_matrix_4.setColumn(0, source_matrix.column(0) * -1);
//    source_matrix_4.setColumn(1, source_matrix.column(1) * -1);
//    source_matrix_4.setColumn(2, source_matrix.column(2));
//    QMatrix4x4 R4 = target_pointcloud.eigvecs_matrix * source_matrix_4.transposed();
//    float trace4 = calculate_angle_from_trace(R4);
//    if (trace4 < trace) {
//        trace = trace4;
//        R = R4;
//    }

    std::vector<QVector3D> new_points;
    QVector3D ct = target_pointcloud.centroid;
    QVector3D cs = source_pointcloud.centroid;

    for (auto point: source_pointcloud.getPoints()) {
        QVector3D new_point = ct + R * (point - cs);
        new_points.push_back(new_point);
    }
    source_pointcloud.set_points(new_points);
    update();
}

float GLWidget::calculate_angle_from_trace(QMatrix4x4 matrix) {
    float trace = matrix.column(0)[0] + matrix.column(1)[1] + matrix.column(2)[2];
    float angle = acos((trace - 1)/2);
    return angle;
}

void GLWidget::radioButton1Clicked()
{
    _projection_algorithm = 0;
    update();
}

void GLWidget::radioButton2Clicked()
{
    _projection_algorithm = 1;
    update();
}

void GLWidget::radioButton3Clicked()
{
    _projection_algorithm = 2;
    update();
}

void GLWidget::checkBox1Clicked() {
    _rays_on = !_rays_on;
    update();
}

void GLWidget::checkBox2Clicked() {
    _projection_on = !_projection_on;
    update();
}

void GLWidget::checkBox3Clicked() {
    _camera1_on = !_camera1_on;
    update();
}

void GLWidget::checkBox4Clicked() {
    _camera2_on = !_camera2_on;
    update();
}

void GLWidget::checkBox5Clicked() {
    _cubes_on = !_cubes_on;
    update();
}

void GLWidget::checkBox6Clicked() {
    draw_kd_tree_partitions = !draw_kd_tree_partitions;
    update();
}

void GLWidget::checkBox7Clicked() {
    draw_oct_tree_partitions = !draw_oct_tree_partitions;
    update();
}

void GLWidget::onDoubleSpinBox1_valueChanged(double arg) {
    _c1_rotation_x = arg;
    update();
}
void GLWidget::onDoubleSpinBox2_valueChanged(double arg) {
    _c1_rotation_y = arg;
    update();
}

void GLWidget::onDoubleSpinBox3_valueChanged(double arg) {
    _c1_rotation_z = arg;
    update();
}

void GLWidget::onDoubleSpinBox4_valueChanged(double arg) {
    _c1_x = arg;
    update();
}

void GLWidget::onDoubleSpinBox5_valueChanged(double arg) {
    _c1_y = arg;
    update();
}

void GLWidget::onDoubleSpinBox6_valueChanged(double arg) {
    _c1_z = arg;
    update();
}

void GLWidget::onDoubleSpinBox7_valueChanged(double arg) {
    _focal_length = arg;
    update();
}

void GLWidget::initShaders()
{
    _shaders.reset(new QOpenGLShaderProgram());
    auto vsLoaded = _shaders->addShaderFromSourceFile(QOpenGLShader::Vertex, ":/vertex_shader.glsl");
    auto fsLoaded = _shaders->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/fragment_shader.glsl");
    assert(vsLoaded && fsLoaded);
    // vector attributes
    _shaders->bindAttributeLocation("vertex", 0);
    _shaders->bindAttributeLocation("pointRowIndex", 1);
    // constants
    _shaders->bind();
    _shaders->setUniformValue("lightPos", QVector3D(0, 0, 50));
    _shaders->setUniformValue("pointsCount", static_cast<GLfloat>(pointcloud.getCount()));
    _shaders->link();
    _shaders->release();

   }

void GLWidget::createContainers()
{
    // create array container and load points into buffer
    const QVector<float>& pointsData =pointcloud.getData();
    //const QVector<float>& pointsData =source_pointcloud.getData();
    //const QVector<float>& pointsData3 =target_pointcloud.getData();
    if(!_vao.isCreated()) _vao.create();
    QOpenGLVertexArrayObject::Binder vaoBinder(&_vao);
    if(!_vertexBuffer.isCreated()) _vertexBuffer.create();
    _vertexBuffer.bind();
    _vertexBuffer.allocate(pointsData.constData(), pointsData.size() * sizeof(GLfloat));
    //_vertexBuffer.allocate(pointsData2.constData(), pointsData2.size() * sizeof(GLfloat));
    //_vertexBuffer.allocate(pointsData3.constData(), pointsData3.size() * sizeof(GLfloat));
    QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
    f->glEnableVertexAttribArray(0);
    f->glEnableVertexAttribArray(1);
    f->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat) + sizeof(GLfloat), nullptr);
    f->glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat) + sizeof(GLfloat), reinterpret_cast<void *>(3*sizeof(GLfloat)));
    _vertexBuffer.release();
}

void GLWidget::setupRenderingCamera()
{
    const RenderingCameraState cameraState = _renderingCamera->state();
    // position and angles
    _cameraMatrix.setToIdentity();
    _cameraMatrix.translate(cameraState.position.x(), cameraState.position.y(), cameraState.position.z());
    _cameraMatrix.rotate   (cameraState.rotation.x(), 1, 0, 0);
    _cameraMatrix.rotate   (cameraState.rotation.y(), 0, 1, 0);
    _cameraMatrix.rotate   (cameraState.rotation.z(), 0, 0, 1);

    // set clipping planes
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    const double rearClippingPlane[] = {0., 0., -1., cameraState.rearClippingDistance};
    glClipPlane(GL_CLIP_PLANE1 , rearClippingPlane);
    const double frontClippingPlane[] = {0., 0., 1., cameraState.frontClippingDistance};
    glClipPlane(GL_CLIP_PLANE2 , frontClippingPlane);

}

void GLWidget::drawPointCloud()
{
    const auto viewMatrix = _projectionMatrix * _cameraMatrix * _worldMatrix;
    _shaders->bind();
    _shaders->setUniformValue("pointsCount", static_cast<GLfloat>(pointcloud.getCount()));
    _shaders->setUniformValue("viewMatrix", viewMatrix);
    _shaders->setUniformValue("pointSize", _pointSize);
    //_shaders->setUniformValue("colorAxisMode", static_cast<GLfloat>(_colorMode));
    _shaders->setUniformValue("colorAxisMode", static_cast<GLfloat>(0));
    _shaders->setUniformValue("pointsBoundMin", pointcloud.getMin());
    _shaders->setUniformValue("pointsBoundMax", pointcloud.getMax());
    glDrawArrays(GL_POINTS, 0, pointcloud.getData().size());
    _shaders->release();
}

void GLWidget::draw_source_pointcloud() {
    draw_points(source_pointcloud.getPoints(), 1, QColor(1.0, 1.0, 1.0));
}

void GLWidget::draw_target_pointcloud() {
    draw_points(target_pointcloud.getPoints(), 1, QColor(1.0, 0.7, 0.0));
}
