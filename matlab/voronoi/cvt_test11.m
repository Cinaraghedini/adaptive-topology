function cvt_test11 ( )

%*****************************************************************************80
%
%% CVT_TEST11 tests CVT.
%
%  Discussion:
%
%    In this test, we initialize the generators to grid points; this is 
%    an unstable CVT solution.  The data would "prefer" to be in a
%    different form.  However, even if we take 2000 steps of CVT iteration,
%    the data is still only slowly progressing towards that other 
%    configuration.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 November 2006
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST11\n' );
  fprintf ( 1, '  CVT computes a Centroidal Voronoi Tessellation.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  In this test, we initialize the generators to\n' );
  fprintf ( 1, '  grid points; this is an unstable CVT solution.\n' );

  dim_num = 2;
  n = 20;
  batch = 1000;
  init = 4;
  init_string = 'user initialization';
  it_max = 60;
  it_fixed = 1;
  sample = 3;
  sample_num = 10000;
  sample_string = 'uniform';
  seed = 123456789;
 
   seed_init = seed;
  
  
  load('Z:\project\database\20\50_16\2\position')

  %load('Z:\Project\code\matlab1608\matlab\NewFailureModel\r')
  
  r=[transpose(position(:,1)); transpose(position(:,2))];
  labels = num2str([1:size(position,1)]');
  xi1=[transpose(r(1,:)) transpose(r(2,:))];
  
  figure(1)
    
  plot(xi1(:,1),xi1(:,2),'LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8,'Marker','o');
  hold all;
  text(xi1(:,1),xi1(:,2), labels, 'VerticalAlignment','bottom','HorizontalAlignment','right','FontName', 'Segoe UI Semibold','Fontsize',15,'FontWeight','bold');

  
  xlim([0 50]);
  ylim([0 50]);

    
  figure(2);
  
  x=r(1,:);
  y=r(2,:);
  
  [vx,vy] = voronoi(x,y);
  
  plot(x,y,'r+',vx,vy,'b-');
  
  xlimit=[min(x)-3 max(x)+3];
  ylimit=([min(y)-3 max(y)+3]);
  
  xlim(xlimit);
  ylim(ylimit);
  
  
  %r8mat_transpose_print ( dim_num, n, r, '  Initial generators (rows):' );

  [ r, seed, it_num, it_diff, energy ] = cvt ( dim_num, n, batch, init, ...
    sample, sample_num, it_max, it_fixed, seed, r );

    figure(3)

    xi=[transpose(r(1,:)) transpose(r(2,:))];
  
    plot(xi(:,1),xi(:,2),'LineStyle','none','MarkerFaceColor',[0.380392163991928 0.380392163991928 0.380392163991928],'MarkerEdgeColor',[0.380392163991928 0.380392163991928 0.380392163991928],'MarkerSize',8,'Marker','o');
    hold all;
    text(xi(:,1),xi(:,2), labels, 'VerticalAlignment','bottom','HorizontalAlignment','right','FontName', 'Segoe UI Semibold','Fontsize',15,'FontWeight','bold');


%   xlim([0 50]);
%   ylim([0 50]);
    
    
    
  figure(4);
  
  x=r(1,:);
  y=r(2,:);
  
  [vx,vy] = voronoi(x,y);
  
  plot(x,y,'r+',vx,vy,'b-');
  
  xlimit=[min(x)-3 max(x)+3];
  ylimit=([min(y)-3 max(y)+3]);
  
  xlim(xlimit);
  ylim(ylimit);
  


  fprintf ( 1, '\n' );
  fprintf ( 1, '  Dimension DIM_NUM =        %12d\n', dim_num );
  fprintf ( 1, '  Number of points N =       %12d\n', n );
  fprintf ( 1, '  Initial SEED =             %12d\n', seed_init );
  fprintf ( 1, '  Current SEED =             %12d\n', seed );
  fprintf ( 1, '  INIT =                    "%s".\n', init_string );
  fprintf ( 1, '  Max iterations IT_MAX =    %12d\n', it_max );
  fprintf ( 1, '  IT_FIXED (fixed samples) = %12d\n', it_fixed );
  fprintf ( 1, '  Iterations IT_NUM =        %12d\n', it_num );
  fprintf ( 1, '  Difference IT_DIFF =       %14f\n', it_diff );
  fprintf ( 1, '  CVT ENERGY =               %14f\n', energy );
  fprintf ( 1, '  SAMPLE =                  "%s".\n', sample_string );
  fprintf ( 1, '  Samples SAMPLE_NUM    =    %12d\n', sample_num );
  fprintf ( 1, '  Sampling BATCH size =      %12d\n', batch );
  fprintf ( 1, '  EPSILON (unit roundoff) =  %12e\n', eps );
  
  r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' );

  return
end