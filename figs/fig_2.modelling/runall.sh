SCRIPT_DIR=~/git/ggr-project/figs/fig_2.modelling
EVAL_DIR="/mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03"
PREDICT_DIR="/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05"
PREDICT_DIR="/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-10-08"

# fig 2b
#$SCRIPT_DIR/fig_2-b.0.plot.encode_clf.R fig_2-b.0 $EVAL_DIR/encode-roadmap.basset.clf.testfold-*/by_task/*txt
#$SCRIPT_DIR/fig_2-b.1.plot.ggr_clf.R fig_2-b.1 $EVAL_DIR/ggr.basset.clf.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.clf.pretrained.folds.testfold-*/by_task/*txt
#$SCRIPT_DIR/fig_2-b.2.plot.ggr_regr.R fig_2-b.1 $EVAL_DIR/ggr.basset.regr.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.regr.pretrained.folds.testfold-*/by_task/*txt

# fig 2c
#$SCRIPT_DIR/fig_2-c.0.plot.scatter.R $EVAL_DIR/ggr.basset.regr.pretrained.folds.testfold-4/by_task/correlation

# fig 2d
#$SCRIPT_DIR/fig_2-d.0.plot.heatmaps.R $PREDICT_DIR/motifs.input_x_grad.early/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.mid/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.late/ggr.scanmotifs.h5
$SCRIPT_DIR/fig_2-d.0.plot.heatmaps.R $PREDICT_DIR/motifs.input_x_grad.dynamic.early/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.dynamic.mid/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.dynamic.late/ggr.scanmotifs.h5
