# 定义函数边界条件
{
  [./topz]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '-1*y+0*z-z+if(t>0,0.5-y,0)' # note that this uses original nodal values of (x,y,z)
  [../]

  [./topz]
    type = FunctionDirichletBC
    variable = z_disp
    boundary = front
    function = 1E-6*t
  [../]

  [./bottomz]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = back
    function = '-1*y+0*z-z+if(t>0,0.5-y,0)' # note that this uses original nodal values of (x,y,z)
  [../]

  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = '5*if(0<t<10,t,1)'
  [../]
}
 
# 输出一条边的位移-时间曲线
