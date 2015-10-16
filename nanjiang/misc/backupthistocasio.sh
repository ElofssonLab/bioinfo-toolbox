#!/bin/bash

rsync -auvz -e ssh --exclude=".*" --exclude="~*" ../bin/ nanjiang@casio.fos.su.se:/data3/wk/MPTopo/bin/
