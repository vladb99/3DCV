#include "mainwindow.h"

#include <QFileDialog>
#include <QKeyEvent>

#include <iostream>
MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent), ui(new Ui::MainWindowClass)
{
    ui->setupUi(this);
	ui->glwidget->setFocusPolicy(Qt::StrongFocus);
    ui->glwidget->setFocus();

    _camera = QSharedPointer<RenderingCamera>(new RenderingCamera());
    ui->glwidget->attachCamera(_camera);

    QObject::connect(ui->pushButton,&QPushButton::clicked,ui->glwidget, &GLWidget::openFileDialog);
    QObject::connect(ui->horizontalSlider, &QSlider::valueChanged, this, &MainWindow::updatePointSize);
    QObject::connect(ui->radioButton_1,&QRadioButton::clicked,ui->glwidget,&GLWidget::radioButton1Clicked);
    QObject::connect(ui->radioButton_2,&QRadioButton::clicked,ui->glwidget,&GLWidget::radioButton2Clicked);
    QObject::connect(ui->radioButton,&QRadioButton::clicked,ui->glwidget,&GLWidget::radioButton3Clicked);
    QObject::connect(ui->checkBox,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox1Clicked);
    QObject::connect(ui->checkBox_2,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox2Clicked);
    QObject::connect(ui->checkBox_3,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox3Clicked);
    QObject::connect(ui->checkBox_4,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox4Clicked);
    QObject::connect(ui->checkBox_5,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox5Clicked);
    QObject::connect(ui->checkBox_6,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox6Clicked);
    QObject::connect(ui->checkBox_7,&QCheckBox::clicked,ui->glwidget,&GLWidget::checkBox7Clicked);
    QObject::connect(ui->pushButton3,&QPushButton::clicked,ui->glwidget, &GLWidget::open_source_pointcloud);
    QObject::connect(ui->pushButton2,&QPushButton::clicked,ui->glwidget, &GLWidget::open_target_pointcloud);
    QObject::connect(ui->pushButton4,&QPushButton::clicked,ui->glwidget, &GLWidget::align_point_clouds);

    updatePointSize(1);

    _camera->setPosition(QVector3D(0.0f, -0.1f, -0.2f));
    _camera->rotate(0, 50, 0);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::updatePointSize(size_t value)
{
    std::cout << "new pointsize: " << value << std::endl;
    ui->glwidget->setPointSize(value);
}

void MainWindow::on_radioButton_clicked()
{

}

void MainWindow::on_doubleSpinBox_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox1_valueChanged(arg1);
}


void MainWindow::on_doubleSpinBox_2_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox2_valueChanged(arg1);
}


void MainWindow::on_doubleSpinBox_3_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox3_valueChanged(arg1);
}


void MainWindow::on_doubleSpinBox_4_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox4_valueChanged(arg1);
}


void MainWindow::on_doubleSpinBox_5_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox5_valueChanged(arg1);
}


void MainWindow::on_doubleSpinBox_6_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox6_valueChanged(arg1);
}


void MainWindow::on_doubleSpinBox_7_valueChanged(double arg1)
{
    ui->glwidget->onDoubleSpinBox7_valueChanged(arg1);
}

void MainWindow::on_checkBox_5_clicked()
{

}

