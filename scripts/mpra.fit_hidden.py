
import os
import h5py
import numpy as np
import pandas as pd

from sklearn import linear_model
from sklearn.model_selection import train_test_split

#from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr

def main():
    """use hidden layer to try and fit to MPRA data
    """
    WORK_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2020-06-29.fit_nn"
    data_file = "{}/ggr.predictions.h5".format(WORK_DIR)

    # pull features
    with h5py.File(data_file, "r") as hf:
        hidden_vals = hf["final_hidden"][:]
        print hidden_vals.shape

    # for each day
    days = ["d0", "d3", "d6"]
    for day in days:
        print day
        with h5py.File(data_file, "r") as hf:
            y = hf[day][:,0]

        # clean a little
        keep_indices = np.where(y != 0)[0]
        X = hidden_vals[keep_indices]
        y = y[keep_indices]

        # splits
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.15, random_state=42)

        print X_train.shape, y_train.shape
        print X_test.shape, y_test.shape
        
        # Create linear regression object
        #regr = linear_model.LinearRegression()
        #regr = linear_model.Ridge(alpha=5)
        #reg = linear_model.RidgeCV(alphas=np.logspace(-6, 6, 13))
        #regr = linear_model.Lasso(alpha=0.005)
        regr = linear_model.LassoCV(cv=5) # <- this is current best
        #regr = linear_model.LassoLarsCV(cv=5)
        #regr = linear_model.ElasticNet(alpha=0.001)
        #regr = linear_model.ElasticNetCV(cv=5, random_state=0)
        #regr = linear_model.Perceptron(tol=1e-3, random_state=0)
        
        # Train the model using the training sets
        regr.fit(X_train, y_train)

        # Make predictions using the testing set
        y_pred = regr.predict(X_test)
        print regr.score(X_test, y_test)

        # save out
        results = pd.DataFrame({
            "MPRA": y_test,
            "NN": y_pred})
            #"NN": y_pred[:,0]})
        results.to_csv("test.txt", sep="\t", header=True, index=False)
        #print pearsonr(results["MPRA"].values, results["NN"].values)
        
        # plot
        script = "/users/dskim89/git/ggr-project/R/plot.mpra.compare_nn.R"
        plot_file = "mpra.compare_nn_hidden.{}.pdf".format(day)
        plot_cmd = "{} {} {}".format(script, "test.txt", plot_file)
        os.system(plot_cmd)
        
    return

main()
