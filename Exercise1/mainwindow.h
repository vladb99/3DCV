//
// (c) Georg Umlauf, 2021
//

#pragma once

#include <QtWidgets/QMainWindow>
#include <QVector3D>
#include <QSharedPointer>

#include "ui_mainwindow.h"

#include "RenderingCamera.h"
class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

protected slots:
  void updatePointSize(size_t);

private slots:
  void on_radioButton_clicked();

  void on_doubleSpinBox_valueChanged(double arg1);

  void on_doubleSpinBox_2_valueChanged(double arg1);

  void on_doubleSpinBox_3_valueChanged(double arg1);

  void on_doubleSpinBox_4_valueChanged(double arg1);

  void on_doubleSpinBox_5_valueChanged(double arg1);

  void on_doubleSpinBox_6_valueChanged(double arg1);

  void on_doubleSpinBox_7_valueChanged(double arg1);

  void on_checkBox_5_clicked();

private:
	Ui::MainWindowClass *ui;
    QSharedPointer<RenderingCamera> _camera;
};
