# Readme

## 文档树介绍

- [Thesis](./ustcthesis/main.pdf)

- `pic` - 所有代码文件
  
  - `warpx` - 利用 [warpx](https://ecp-warpx.github.io/) 的模拟代码 （主要代码内容）
    
    - `run_3d` - 三维模拟
    
    - `run_2d` - 二维模拟
    
    - `script` - 超算相关运行文件
    
    - `test` - 测试相关代码
  
  - `smilei` - 利用 [smilei](https://smileipic.github.io/Smilei/) 的模拟代码
  
  - `Archive` - 主要包含的是测试代码即输入文件, 删除了所有输出文件(文件夹名大概反映了相关的代码内容, 如 `injection` 表示测试的是和粒子注入相关的内容)
  
  - `plot` - 画图相关的部分方便调用的代码

## 其他说明

模拟版本采用了简单的数字来区分, 如 `01` 表示第一次成功模拟的代码, 不同版本的代码区别在于参数的调整/代码的的重构, 建议从最新的代码版本开始看起 (如`run2d`中的`06`), 利用 `diff` 或者代码编辑器查看不同版本的区别.

所有输出文件已被删除(每个输出文件占用空间在GB量级), 可利用 `.pbs` 执行相关 `.py` 或 `input` 文件重新获得输出.

每个任务文件夹中包含的 `post_processing.*` 文件表示的是画图或分析相关的代码; 三维可视化采用的是 `paraview` 手工逐步调整未保存相关代码.

#### 联系方式

张子衿, zijin@ucla.edu


