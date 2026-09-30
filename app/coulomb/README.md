# 库仑特殊函数 COUL90

COUL90 是 A. R. Barnett 的双精度 Fortran 程序。本目录保留来源记录与固定校验值；公开仓库不包含 `coul90.f`。来源为 [FRESCO 作者维护页面](https://www.fresco.org.uk/functions/barnett/index.htm)，具体下载地址记录在 `UPSTREAM.json` 的 `source_url`。

## 安装与校验

在程序包根目录运行 `make`、`make verify-sources` 或 `make check` 时，缺失的文件由 `fetch-external` 目标自动获取。也可先单独运行：

```sh
make fetch-external
```

安装脚本需要 Bash、GNU 基础命令以及 curl 或 wget。HTTPS 下载保持证书验证。下载先写入临时文件，SHA-256 符合下列固定值后才安装为 `app/coulomb/coul90.f`；已有匹配文件直接通过。下载或校验失败时保留原文件并返回失败。

```text
50ad795901114080ded8f63ebe3895b8460865f10fd3e59ff39d60e48456d13f
```

离线环境可从原站取得同一版本文件，使用：

```sh
./fetch_external.sh --source /path/to/COUL90.FOR
```

离线拷入采用同一校验值。文件按原始字节保存，包括原有换行；接口、算法和参数均沿用已验证的源文件。文件名已加入 `.gitignore`，本地下载不纳入公开仓库。

## 调用与来源范围

`pp_scattering.f90` 仅调用 `KFN=0`，获得正能量实库仑函数 F_L、G_L 及对 x 的导数；eta=0 时给出自由 Riccati 函数。外层封装检查 IFAIL、有限值和 Wronskian。COUL90 的 COMMON 区使此例程不能同时由多个线程调用；本项目散射程序按序调用。

编译时按固定格式 Fortran 文件加入命令行，不依赖 Python、C 或额外动态特殊函数库。记录的原来源没有明确许可证文本，本项目不为 COUL90 另行授予许可。软件根许可证不覆盖这份由用户从原站获取的第三方源文件；使用与后续分发应依据其来源条件。
