function centeredattack_img=crop(img,percentage)
    % 裁剪攻击    
    
        img_cc=img;
        percentage = 0.01 * percentage;
        [row,col]=size(img_cc);% row为图像的长,col为图像的宽
        row_begin=1;
        row_end=round(row*sqrt(percentage));
        col_begin = 1;
        col_end = round(col*sqrt(percentage));
        img_cc(row_begin:row_end,col_begin:col_end,:) = 0;% 将选中的区域的值置为0,即显示为黑色
        centeredattack_img=img_cc;
    end