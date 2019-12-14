function [cartesian_cell] = cell_cartesian(cell_1,cell_2)

%takes cartesian product of two cells and puts them in a grid.  There has
%to be a better way to do this

[dim_1,dim_2] = meshgrid(1:length(cell_1),1:length(cell_2))
dim_1 = num2cell(dim_1');
dim_2 = num2cell(dim_2');

%cartesian_cell = repmat({{'',''}},length(cell_1),length(cell_2))

cartesian_cell = cellfun(@(x,y) {cell_1{x},cell_2{y}} ,dim_1,dim_2,'UniformOutput',false)



end

