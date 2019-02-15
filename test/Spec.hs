{-# language TemplateHaskell #-}

module Main (main) where

import Data.Complex (Complex(..))

import Data.Primitive.PrimArray (PrimArray)
import qualified Data.Primitive.PrimArray as PA

import Hedgehog
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range

import Data.Primitive.Instances ()

import qualified Data.Primitive.Contiguous.FFT as CF
import qualified Numeric.MathFunctions.Comparison as MF

main :: IO Bool
main = checkParallel $$(discover)

prop_close :: Property
prop_close = property $ do
  x <- forAll genPrimArray
  let fft' = fft x
  let ifft' = ifft fft'
  pwithin 5 ifft' x === True

fft,ifft :: PrimArray (Complex Double) -> PrimArray (Complex Double)
fft = CF.fft
ifft = CF.ifft

pwithin :: Int -> PrimArray (Complex Double) -> PrimArray (Complex Double) -> Bool
pwithin acc xs ys = all (==True) (zipWith (within acc) (PA.primArrayToList xs) (PA.primArrayToList ys))

within :: Int -> Complex Double -> Complex Double -> Bool
within acc (x :+ y) (x' :+ y') = MF.within acc x x' && MF.within acc y y'

genDouble :: Gen Double
genDouble = Gen.double (Range.linearFrac 5 200)

genComplex :: Gen (Complex Double)
genComplex = (:+) <$> genDouble <*> genDouble

genPrimArray :: Gen (PrimArray (Complex Double))
genPrimArray = PA.primArrayFromList <$> Gen.list (Range.singleton 8) genComplex
