![](Aspose.Words.763d9e3f-0fe8-4e4f-8ab0-f990fed13797.001.png)

静磁場方程式 

Ñ×***B***=0, ***B***=*m**H***, Ñ´***H*** =0 (1) 

磁気スカラポテンシャル

***H*** =Ñ*w* (2) 支配方程式 

Ñ×*m*Ñ*w* =0 *in* W (3) 境界条件 

- =*ws on* G *D*

*m*Ñ*w*×***n***= ***B***×***n***= *Bs on* G*N* (4) G*D* +G*N* =¶W

重み付き残差法 

*W* = -ò *y*Ñ×*m*Ñ*w* = ò Ñ*y* ×*m*Ñ*w* -ò Ñ×(*ym*Ñ*w*)

- ò Ñ*y*W ×*m*Ñ*w* -ò ***n***W×(*ym*Ñ*w*) W
  - ¶W
- ò Ñ*y* ×*m*Ñ*w* -ò *y Bs* (5) 

W G*N*

- = 0 *on* G*D*

離散化 

- å å
- = *wnNn* + *wnNn* + *wsnNn*

*n*ÎW *n*ÎG*N n*ÎG*D*

- = *Nm*, *m*ÎW+G*N*

*W* = ò Ñ*y* ×*m*Ñ*w* -ò *yBs*

åW G*N* å

- *wn* ò Ñ*Nm* ×*m*Ñ*Nn* + *wsn* ò Ñ*Nm* ×*m*Ñ*Nn* -ò *NmBs n*ÎW+G*N* W *n*ÎG*D* W G*N*
- å

*wn* ò Ñ*Nm* ×*m*Ñ*Nn* = ò *NmBs* - *wsn* ò Ñ*Nm* ×*m*Ñ*Nn n*ÎW+G*N* W G*N n*ÎG*D* W
