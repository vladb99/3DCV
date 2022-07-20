//
// Widget für Interaktion und Kontrolle
//
// (c) Georg Umlauf, 2021
//

#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QMatrix4x4>
#include <QMatrix3x3>
#include <QVector3D>
#include <QSharedPointer>

#include <vector>

#include "RenderingCamera.h"
#include "PointCloud.h"
#include "camera.h"
#include "cube.h"
#include "kdtree.h"
#include "kdnode.h"
#include "octtree.h"
#include "octnode.h"
#include <math.h>

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    GLWidget(QWidget* parent = nullptr);
    ~GLWidget() Q_DECL_OVERRIDE;

public slots:
    // open a PLY file
    void openFileDialog();
    void radioButton1Clicked();
    void radioButton2Clicked();
    void radioButton3Clicked();
    void checkBox1Clicked();
    void checkBox2Clicked();
    void checkBox3Clicked();
    void checkBox4Clicked();
    void checkBox5Clicked();    
    void checkBox6Clicked();
    void checkBox7Clicked();
    void onDoubleSpinBox1_valueChanged(double arg);
    void onDoubleSpinBox2_valueChanged(double arg);
    void onDoubleSpinBox3_valueChanged(double arg);
    void onDoubleSpinBox4_valueChanged(double arg);
    void onDoubleSpinBox5_valueChanged(double arg);
    void onDoubleSpinBox6_valueChanged(double arg);
    void onDoubleSpinBox7_valueChanged(double arg);
    void setPointSize(size_t size);
    void attachCamera(QSharedPointer<RenderingCamera> camera);
    void open_source_pointcloud();
    void open_target_pointcloud();
    void align_point_clouds();

protected:
    void paintGL() Q_DECL_OVERRIDE;
    void initializeGL() Q_DECL_OVERRIDE;
    void resizeGL(int width, int height) Q_DECL_OVERRIDE;

    // navigation
    void keyPressEvent(QKeyEvent   *event) Q_DECL_OVERRIDE;
    void wheelEvent(QWheelEvent *) Q_DECL_OVERRIDE;
    void mouseMoveEvent (QMouseEvent *event) Q_DECL_OVERRIDE;


private slots:
    void onCameraChanged(const RenderingCameraState& state);

private:
    int global_index = 0;
    // interaction control
    bool X_Pressed=false, Y_Pressed=false;
    QPoint _prevMousePosition;

    bool test = true;

    // shader control
    void initShaders();
    void createContainers();
    QOpenGLVertexArrayObject _vao;
    QOpenGLBuffer _vertexBuffer;
    QScopedPointer<QOpenGLShaderProgram> _shaders;

    // rendering control
    void setupRenderingCamera();
    QSharedPointer<RenderingCamera> _renderingCamera;
    QMatrix4x4 _projectionMatrix;
    QMatrix4x4 _cameraMatrix;
    QMatrix4x4 _worldMatrix;

    std::vector<std::pair<QVector3D, QColor> > _axesLines;

    // algorithm control
    int _projection_algorithm = 0;

    // projection control
    bool _rays_on = false;
    bool _projection_on = true;
    bool _camera1_on = true;
    bool _camera2_on = true;
    bool _cubes_on = true;

    // camera1 control
    float _c1_rotation_x = 0;
    float _c1_rotation_y = 0;
    float _c1_rotation_z = 0;

    float _c1_x = 4;
    float _c1_y = 0;
    float _c1_z = 0;

    float _focal_length = 2;

    // partitions
    bool draw_kd_tree_partitions = false;
    bool draw_oct_tree_partitions = false;

    PointCloud pointcloud;
    KdTree kdtree;
    OctTree octtree;
    OctNode *octtree_root;;

    // scene and scene control
    void cleanup       ();
    void drawScene     ();
    void drawPointCloud();
    void drawFrameAxis ();
    float _pointSize;

    // PCA
    PointCloud source_pointcloud;
    PointCloud target_pointcloud;

    QVector4D rotate_point_y_axis(QVector4D vector, float angle);
    QVector4D rotate_point_x_axis(QVector4D vector, float angle);
    QVector4D rotate_point_z_axis(QVector4D vector, float angle);

    std::vector<QVector3D> get_cube_points(QVector3D start_position, float size, float angleX, float angleY, float angleZ ,QColor color);

    void draw_cube(Cube cube, QColor color);
    void draw_cube_2(Cube cube, QColor color);
    void draw_camera_model(Camera camera, float height, float width);
    void draw_projection(std::vector<QVector3D> projected_points);
    void draw_rays(Camera camera, Cube cube);
    void perspective_reconstruction(Camera left_camera, Camera right_camera, Cube cube);

    void construct_kd_tree(int l_index, int r_index);
    void draw_plane(QVector3D point, int direction, double size);
    void draw_points(std::vector<QVector3D> points, int size, QColor color);

    void draw_pca_axes(PointCloud* pc);
    float multiply_rows(std::vector<float> row1, std::vector<float> row2);
    std::vector<float> eigvals(QMatrix4x4 matrix);
    std::vector<QVector3D> eigvecs(QMatrix4x4 matrix, std::vector<float> vals);
    std::vector<float> solve_cubic_equation(float a, float b, float c, float d);

    void draw_axes(std::vector<QVector3D> axes, QVector3D centroid, float scale);
    void draw_point(QVector3D point);

    void draw_source_pointcloud();
    void draw_target_pointcloud();

    float calculate_angle_from_trace(QMatrix4x4 matrix);
    void update_source_pointcloud();
};
