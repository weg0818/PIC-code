
## gbc-master

## 恢复libs中的子库
1. git submodule update --init --recursive
--init：初始化 .gitmodules 中列出的子模块。
--recursive：如果子模块内部还嵌套了子模块（递归），也会一并下载。

2. 调整位置
如 gbc-master/2d/extra/libs/external/tetgen 移动到  gbc-master/2d/extra/libs/libigl/external/tetgen


