function QI_matrix = QI2matrix_HS(QI_, Method_, Time_)
k = length(QI_);

col = 1;
% First row
QI_matrix{1,col} = 'Method'; col = col+1;
QI_matrix{1,col} = 'CC';col = col+1;
QI_matrix{1,col} = 'SAM';col = col+1;
QI_matrix{1,col} = 'RMSE';col = col+1;
QI_matrix{1,col} = 'ERGAS';col = col+1;
QI_matrix{1,col} = 'PSNR';col = col+1;
QI_matrix{1,col} = 'Time';col = col+1;

for i = 1:k
    col = 1;
    QI_matrix{i+1,col} = Method_{i}; col = col+1;
    QI_matrix{i+1,col} = QI_{i}.cc; col = col+1;
    QI_matrix{i+1,col} = QI_{i}.sam; col = col+1;
    QI_matrix{i+1,col} = QI_{i}.rmse; col = col+1;
    QI_matrix{i+1,col} = QI_{i}.ergas; col = col+1;
    QI_matrix{i+1,col} = QI_{i}.psnr; col = col+1;    
    QI_matrix{i+1,col} = Time_{i}; col = col+1;
end