class ArrayFunctionTransformer(sk.base.BaseEstimator, sk.base.TransformerMixin):
    def __init__(self, FUN, target, params=None):
        self.FUN = FUN
        self.target = target
        self.params = params
    def fit(self,X,y=None):
        return self
    def transform(self,X):
        X = np.asarray(getattr(X,self.target))
        result = []
        for x in X:
            result.append(self.FUN(x, *self.params if self.params else ()))
        return np.asarray(result)