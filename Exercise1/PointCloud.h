#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <QString>
#include <QVector>
#include <QVector3D>
#include <vector>
#include <QMatrix4x4>

static const size_t POINT_STRIDE = 4; // x, y, z, index

class PointCloud
{
public:
    PointCloud();
    ~PointCloud();

    bool loadPLY(const QString&);

private:
    size_t _pointsCount = 0;
    QVector<float> _pointsData;
    QVector3D _pointsBoundMin;
    QVector3D _pointsBoundMax;
    std::vector<QVector3D> _points;

public:
    size_t getCount() const { return _pointsCount; }
    QVector3D getMin() const { return _pointsBoundMin; }
    QVector3D getMax() const { return _pointsBoundMax; }
    const QVector<float>& getData() const { return _pointsData; }
    std::vector<QVector3D> getPoints() const { return _points; }
    void set_points(std::vector<QVector3D> points) { _points = points; }

    // sorted points
    std::vector<QVector3D> _X;
    std::vector<QVector3D> _Y;
    std::vector<QVector3D> _Z;

    std::vector<QVector3D> _X_q;
    std::vector<QVector3D> _Y_q;
    std::vector<QVector3D> _Z_q;

    // PCA
    QVector3D centroid;
    QMatrix4x4 eigvecs_matrix;

    bool test = false;

    struct sort_after_x
    {
        inline bool operator() (const QVector3D point1, const QVector3D point2)
        {
            return (point1.x() < point2.x());
        }
    };

    struct sort_after_y
    {
        inline bool operator() (const QVector3D point1, const QVector3D point2)
        {
            return (point1.y() < point2.y());
        }
    };

    struct sort_after_z
    {
        inline bool operator() (const QVector3D point1, const QVector3D point2)
        {
            return (point1.z() < point2.z());
        }
    };
};

#endif // POINTCLOUD_H
