USE teglon;

-- MySQL dump 10.13  Distrib 8.0.12, for macos10.13 (x86_64)
--
-- Host: 127.0.0.1    Database: test_teglon
-- ------------------------------------------------------
-- Server version	8.0.16

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
 SET NAMES utf8 ;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `Band`
--

DROP TABLE IF EXISTS `Band`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `Band` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(45) NOT NULL,
  `Effective_Wavelength` double NOT NULL,
  `F99_Coefficient` double NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=14 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Detector`
--

DROP TABLE IF EXISTS `Detector`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `Detector` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(45) NOT NULL,
  `Deg_width` double DEFAULT NULL,
  `Deg_height` double DEFAULT NULL,
  `Deg_radius` double DEFAULT NULL,
  `Area` double DEFAULT NULL,
  `MinDec` double NOT NULL DEFAULT '-90',
  `MaxDec` double NOT NULL DEFAULT '90',
  `Poly` multipolygon DEFAULT NULL,
  `TM_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=54 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Distance`
--

DROP TABLE IF EXISTS `Distance`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `Distance` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Galaxy_id` int(11) NOT NULL,
  `z` double NOT NULL,
  `z_dist` double NOT NULL,
  `z_dist_err` double NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Distance_Galaxy_id_idx` (`Galaxy_id`),
  KEY `idx_Distance_z` (`z`),
  CONSTRAINT `Distance_Galaxy_id` FOREIGN KEY (`Galaxy_id`) REFERENCES `Galaxy` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=2721740 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GWGC`
--

DROP TABLE IF EXISTS `GWGC`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `GWGC` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `PGC` int(11) DEFAULT NULL,
  `Name` varchar(64) DEFAULT NULL,
  `RA` double NOT NULL,
  `Dec` double NOT NULL,
  `TT` double DEFAULT NULL,
  `B` double DEFAULT NULL,
  `e_B` double DEFAULT NULL,
  `semi_a` double DEFAULT NULL,
  `e_semi_a` double DEFAULT NULL,
  `semi_b` double DEFAULT NULL,
  `e_semi_b` double DEFAULT NULL,
  `b_a` double DEFAULT NULL,
  `e_b_a` double DEFAULT NULL,
  `PA` double DEFAULT NULL,
  `B_abs` double DEFAULT NULL,
  `e_B_abs` double DEFAULT NULL,
  `Dist` double DEFAULT NULL,
  `e_Dist` double DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=106625 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Galaxy`
--

DROP TABLE IF EXISTS `Galaxy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `Galaxy` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `PGC` int(11) DEFAULT NULL,
  `Name_GWGC` varchar(64) DEFAULT NULL,
  `Name_HyperLEDA` varchar(64) DEFAULT NULL,
  `Name_2MASS` varchar(64) DEFAULT NULL,
  `Name_SDSS_DR12` varchar(45) DEFAULT NULL,
  `RA` double NOT NULL,
  `_Dec` double NOT NULL,
  `dist` double DEFAULT NULL,
  `dist_err` double DEFAULT NULL,
  `z` double DEFAULT NULL,
  `B` double DEFAULT NULL,
  `B_err` double DEFAULT NULL,
  `B_abs` double DEFAULT NULL,
  `J` double DEFAULT NULL,
  `J_err` double DEFAULT NULL,
  `H` double DEFAULT NULL,
  `H_err` double DEFAULT NULL,
  `K` double DEFAULT NULL,
  `K_err` double DEFAULT NULL,
  `flag1` varchar(8) DEFAULT NULL,
  `flag2` varchar(8) DEFAULT NULL,
  `flag3` varchar(8) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `idx_Galaxy_RA` (`RA`),
  KEY `idx_Galaxy__Dec` (`_Dec`),
  KEY `idx_Galaxy_B` (`B`),
  KEY `idx_Galaxy_B_abs` (`B_abs`),
  KEY `idx_Galaxy_dist` (`dist`)
) ENGINE=InnoDB AUTO_INCREMENT=3262883 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GalaxyDistance`
--

DROP TABLE IF EXISTS `GalaxyDistance`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `GalaxyDistance` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Galaxy_id` int(11) NOT NULL,
  `Distance_id` int(11) NOT NULL,
  `PGC` int(11) DEFAULT NULL,
  `Name_GWGC` varchar(64) DEFAULT NULL,
  `Name_HyperLEDA` varchar(64) DEFAULT NULL,
  `Name_2MASS` varchar(64) DEFAULT NULL,
  `Name_SDSS_DR12` varchar(45) DEFAULT NULL,
  `RA` double NOT NULL,
  `_Dec` double NOT NULL,
  `dist` double DEFAULT NULL,
  `dist_err` double DEFAULT NULL,
  `z_dist` double DEFAULT NULL,
  `z_dist_err` double DEFAULT NULL,
  `z` double DEFAULT NULL,
  `B` double DEFAULT NULL,
  `B_err` double DEFAULT NULL,
  `B_abs` double DEFAULT NULL,
  `J` double DEFAULT NULL,
  `J_err` double DEFAULT NULL,
  `H` double DEFAULT NULL,
  `H_err` double DEFAULT NULL,
  `K` double DEFAULT NULL,
  `K_err` double DEFAULT NULL,
  `flag1` varchar(8) DEFAULT NULL,
  `flag2` varchar(8) DEFAULT NULL,
  `flag3` varchar(8) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `idx_GalaxyDistance_id` (`id`),
  KEY `idx_GalaxyDistance_Galaxy_id` (`Galaxy_id`),
  KEY `idx_GalaxyDistance_Distance_id` (`Distance_id`),
  KEY `idx_GalaxyDistance_RA` (`RA`),
  KEY `idx_GalaxyDistance__Dec` (`_Dec`),
  KEY `idx_GalaxyDistance_dist` (`dist`),
  KEY `idx_GalaxyDistance_z_dist` (`z_dist`),
  KEY `idx_GalaxyDistance_z_dist_err` (`z_dist_err`),
  KEY `idx_GalaxyDistance_z` (`z`),
  KEY `idx_GalaxyDistance_B` (`B`)
) ENGINE=InnoDB AUTO_INCREMENT=2752471 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GalaxyDistance2`
--

DROP TABLE IF EXISTS `GalaxyDistance2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `GalaxyDistance2` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Galaxy_id` int(11) NULL,
  `Distance_id` int(11) NULL,
  `PGC` int(11) DEFAULT NULL,
  `Name_GWGC` varchar(64) DEFAULT NULL,
  `Name_HyperLEDA` varchar(64) DEFAULT NULL,
  `Name_2MASS` varchar(64) DEFAULT NULL,
  `Name_SDSS_DR12` varchar(64) DEFAULT NULL,
  `RA` double NOT NULL,
  `_Dec` double NOT NULL,
  `Coord` point NOT NULL,
  `dist` double DEFAULT NULL,
  `dist_err` double DEFAULT NULL,
  `z_dist` double DEFAULT NULL,
  `z_dist_err` double DEFAULT NULL,
  `z` double DEFAULT NULL,
  `B` double DEFAULT NULL,
  `B_err` double DEFAULT NULL,
  `B_abs` double DEFAULT NULL,
  `J` double DEFAULT NULL,
  `J_err` double DEFAULT NULL,
  `H` double DEFAULT NULL,
  `H_err` double DEFAULT NULL,
  `K` double DEFAULT NULL,
  `K_err` double DEFAULT NULL,
  `flag1` varchar(8) DEFAULT NULL,
  `flag2` varchar(8) DEFAULT NULL,
  `flag3` varchar(8) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `idx_GalaxyDistance2_id` (`id`),
  KEY `idx_GalaxyDistance2_Galaxy_id` (`Galaxy_id`),
  KEY `idx_GalaxyDistance2_Distance_id` (`Distance_id`),
  KEY `idx_GalaxyDistance2_RA` (`RA`),
  KEY `idx_GalaxyDistance2__Dec` (`_Dec`),
  KEY `idx_GalaxyDistance2_dist` (`dist`),
  KEY `idx_GalaxyDistance2_z_dist` (`z_dist`),
  KEY `idx_GalaxyDistance2_z_dist_err` (`z_dist_err`),
  KEY `idx_GalaxyDistance2_z` (`z`),
  KEY `idx_GalaxyDistance2_B` (`B`),
  SPATIAL KEY `idx_GalaxyDistance2_Coord` (`Coord`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `HealpixMap`
--

DROP TABLE IF EXISTS `HealpixMap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `HealpixMap` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `GWID` varchar(45) NOT NULL,
  `URL` text NOT NULL,
  `Filename` varchar(255) NOT NULL,
  `NSIDE` int(11) NOT NULL,
  `t_0` double NOT NULL,
  `SubmissionTime` datetime NOT NULL,
  `NetProbToGalaxies` double NOT NULL,
  `RescaledNSIDE` int(11) NOT NULL DEFAULT '128',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=10 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `HealpixPixel`
--

DROP TABLE IF EXISTS `HealpixPixel`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `HealpixPixel` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `HealpixMap_id` int(11) NOT NULL,
  `Pixel_Index` int(11) NOT NULL,
  `Prob` double NOT NULL,
  `Distmu` double NOT NULL,
  `Distsigma` double NOT NULL,
  `Distnorm` double NOT NULL,
  `Mean` double NOT NULL,
  `Stddev` double NOT NULL,
  `Norm` double NOT NULL,
  `N128_SkyPixel_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `idx_HealpixPixel_N32_SkyPixel_id` (`N128_SkyPixel_id`),
  KEY `idx_HealpixPixel_Prob` (`Prob`),
  KEY `idx_HealpixPixel_Pixel_Index` (`Pixel_Index`),
  KEY `FK_HealpixPixel_HealpixMap_idx` (`HealpixMap_id`),
  CONSTRAINT `FK_HealpixPixel_HealpixMap` FOREIGN KEY (`HealpixMap_id`) REFERENCES `HealpixMap` (`id`),
  CONSTRAINT `FK_HealpixPixel_SkyPixel` FOREIGN KEY (`N128_SkyPixel_id`) REFERENCES `SkyPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=68222824 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `HealpixPixel_Completeness`
--

DROP TABLE IF EXISTS `HealpixPixel_Completeness`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `HealpixPixel_Completeness` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `HealpixPixel_id` int(11) NOT NULL,
  `PixelCompleteness` double NOT NULL,
  `Renorm2DProb` double NOT NULL,
  `NetPixelProb` double NOT NULL,
  `HealpixMap_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_HealpixPixel_Completeness_HealpixPixel_idx` (`HealpixPixel_id`),
  KEY `FK_HealpixPixel_Completeness_HealpixMap_id_idx` (`HealpixMap_id`),
  KEY `idx_HealpixPixel_Completeness_NetPixelProb` (`NetPixelProb`),
  CONSTRAINT `FK_HealpixPixel_Completeness_HealpixMap_id` FOREIGN KEY (`HealpixMap_id`) REFERENCES `HealpixMap` (`id`),
  CONSTRAINT `FK_HealpixPixel_Completeness_HealpixPixel` FOREIGN KEY (`HealpixPixel_id`) REFERENCES `HealpixPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=64225141 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `HealpixPixel_GalaxyDistance2`
--

DROP TABLE IF EXISTS `HealpixPixel_GalaxyDistance2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `HealpixPixel_GalaxyDistance2` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `HealpixPixel_id` int(11) NOT NULL,
  `GalaxyDistance2_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_HealpixPixel_GalaxyDistance2_HealpixPixel_id_idx` (`HealpixPixel_id`),
  KEY `FK_HealpixPixel_GalaxyDistance2_GalaxyDistance2_id_idx` (`GalaxyDistance2_id`),
  KEY `idx_HealpixPixel_GalaxyDistance2_GalaxyDistance2_id` (`GalaxyDistance2_id`),
  CONSTRAINT `FK_HealpixPixel_GalaxyDistance2_GalaxyDistance2_id` FOREIGN KEY (`GalaxyDistance2_id`) REFERENCES `GalaxyDistance2` (`id`),
  CONSTRAINT `FK_HealpixPixel_GalaxyDistance2_HealpixPixel_id` FOREIGN KEY (`HealpixPixel_id`) REFERENCES `HealpixPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=23282886 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `HealpixPixel_GalaxyDistance2_Weight`
--

DROP TABLE IF EXISTS `HealpixPixel_GalaxyDistance2_Weight`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `HealpixPixel_GalaxyDistance2_Weight` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `HealpixPixel_GalaxyDistance2_id` int(11) NOT NULL,
  `LumWeight` double NOT NULL,
  `zWeight` double NOT NULL,
  `Prob2DWeight` double NOT NULL,
  `Norm4DWeight` double NOT NULL,
  `GalaxyProb` double NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_HGW_HealpixPixel_GalaxyDistance2_id_idx` (`HealpixPixel_GalaxyDistance2_id`),
  KEY `idx_HealpixPixel_GalaxyDistance2_Weight_GalaxyProb` (`GalaxyProb`),
  CONSTRAINT `FK_HGW_HealpixPixel_GalaxyDistance2_id` FOREIGN KEY (`HealpixPixel_GalaxyDistance2_id`) REFERENCES `HealpixPixel_GalaxyDistance2` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=23282886 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ObservedTile`
--

DROP TABLE IF EXISTS `ObservedTile`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `ObservedTile` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `FileName` varchar(256) DEFAULT NULL,
  `Detector_id` int(11) NOT NULL,
  `FieldName` varchar(45) NOT NULL,
  `RA` double NOT NULL,
  `_Dec` double NOT NULL,
  `Coord` point NOT NULL,
  `Poly` multipolygon NOT NULL,
  `EBV` double NOT NULL,
  `N128_SkyPixel_id` int(11) NOT NULL,
  `Band_id` int(11) NOT NULL,
  `MJD` double DEFAULT NULL,
  `Exp_Time` double DEFAULT NULL,
  `Mag_Lim` double DEFAULT NULL,
  `HealpixMap_id` int(11) NOT NULL,
  `PositionAngle` double NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `FK_ObservedTile_Band_Band_id_idx` (`Band_id`,`HealpixMap_id`,`Detector_id`,`N128_SkyPixel_id`),
  KEY `FK_ObservedTile_Detector_Detector_id` (`Detector_id`),
  KEY `FK_ObservedTile_Detector_SkyPixel_id` (`N128_SkyPixel_id`),
  KEY `FK_ObservedTile_Detector_HealpixMap_id` (`HealpixMap_id`),
  CONSTRAINT `FK_ObservedTile_Detector_Band_id` FOREIGN KEY (`Band_id`) REFERENCES `Band` (`id`),
  CONSTRAINT `FK_ObservedTile_Detector_Detector_id` FOREIGN KEY (`Detector_id`) REFERENCES `Detector` (`id`),
  CONSTRAINT `FK_ObservedTile_Detector_HealpixMap_id` FOREIGN KEY (`HealpixMap_id`) REFERENCES `HealpixMap` (`id`),
  CONSTRAINT `FK_ObservedTile_Detector_SkyPixel_id` FOREIGN KEY (`N128_SkyPixel_id`) REFERENCES `SkyPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=21885 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ObservedTile_HealpixPixel`
--

DROP TABLE IF EXISTS `ObservedTile_HealpixPixel`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `ObservedTile_HealpixPixel` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `ObservedTile_id` int(11) NOT NULL,
  `HealpixPixel_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_ObservedTile_HealpixPixel_ObservedTile_id_idx` (`ObservedTile_id`),
  KEY `FK_ObservedTile_HealpixPixel_HealpixPixel_id_idx` (`HealpixPixel_id`),
  CONSTRAINT `FK_ObservedTile_HealpixPixel_HealpixPixel_id` FOREIGN KEY (`HealpixPixel_id`) REFERENCES `HealpixPixel` (`id`),
  CONSTRAINT `FK_ObservedTile_HealpixPixel_ObservedTile_id` FOREIGN KEY (`ObservedTile_id`) REFERENCES `ObservedTile` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=950769 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `SkyCompleteness`
--

DROP TABLE IF EXISTS `SkyCompleteness`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `SkyCompleteness` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `SkyPixel_id` int(11) NOT NULL,
  `SkyDistance_id` int(11) NOT NULL,
  `L10` double NOT NULL,
  `Completeness` double NOT NULL,
  `SmoothedCompleteness` double NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `Completeness_SkyPixel_id_idx` (`SkyPixel_id`),
  KEY `Completeness_SkyDistance_id_idx` (`SkyDistance_id`),
  CONSTRAINT `Completeness_SkyDistance_id` FOREIGN KEY (`SkyDistance_id`) REFERENCES `SkyDistance` (`id`),
  CONSTRAINT `Completeness_SkyPixel_id` FOREIGN KEY (`SkyPixel_id`) REFERENCES `SkyPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=2549809 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `SkyDistance`
--

DROP TABLE IF EXISTS `SkyDistance`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `SkyDistance` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `D1` double NOT NULL,
  `D2` double NOT NULL,
  `dCoV` double NOT NULL,
  `Sch_L10` double NOT NULL,
  `NSIDE` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `idx_SkyDistance_D1` (`D1`),
  KEY `idx_SkyDistance_D2` (`D2`)
) ENGINE=InnoDB AUTO_INCREMENT=277 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `SkyPixel`
--

DROP TABLE IF EXISTS `SkyPixel`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `SkyPixel` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `RA` double NOT NULL,
  `_Dec` double NOT NULL,
  `Coord` point NOT NULL,
  `Poly` multipolygon NOT NULL,
  `NSIDE` int(11) NOT NULL,
  `Pixel_Index` int(11) NOT NULL,
  `Parent_Pixel_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `idx_SkyPixel_RA` (`RA`),
  KEY `idx_SkyPixel__Dec` (`_Dec`),
  SPATIAL KEY `idx_SkyPixel_Coord` (`Coord`),
  SPATIAL KEY `idx_SkyPixel_Poly` (`Poly`),
  KEY `idx_SkyPixel_Pixel_Index` (`Pixel_Index`),
  KEY `idx_SkyPixel_NSIDE` (`NSIDE`),
  KEY `FK_SkyPixel_SkyPixel_Parent1` (`Parent_Pixel_id`),
  CONSTRAINT `FK_SkyPixel_SkyPixel_Parent1` FOREIGN KEY (`Parent_Pixel_id`) REFERENCES `SkyPixel` (`id`) ON DELETE RESTRICT ON UPDATE RESTRICT
) ENGINE=InnoDB AUTO_INCREMENT=262129 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `SkyPixel_EBV`
--

DROP TABLE IF EXISTS `SkyPixel_EBV`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `SkyPixel_EBV` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `N128_SkyPixel_id` int(11) NOT NULL,
  `EBV` double NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_SkyPixel_EBV_N128_SkyPixel_id_idx` (`N128_SkyPixel_id`),
  CONSTRAINT `FK_SkyPixel_EBV_N128_SkyPixel_id` FOREIGN KEY (`N128_SkyPixel_id`) REFERENCES `SkyPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=196609 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `SkyPixel_GalaxyDistance2`
--

DROP TABLE IF EXISTS `SkyPixel_GalaxyDistance2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `SkyPixel_GalaxyDistance2` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `SkyPixel_id` int(11) NOT NULL,
  `GalaxyDistance2_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_GalaxyDistance2_id_idx` (`GalaxyDistance2_id`),
  KEY `FK_SkyPixel_id_idx` (`SkyPixel_id`),
  CONSTRAINT `FK_GalaxyDistance2_id` FOREIGN KEY (`GalaxyDistance2_id`) REFERENCES `GalaxyDistance2` (`id`),
  CONSTRAINT `FK_SkyPixel_id` FOREIGN KEY (`SkyPixel_id`) REFERENCES `SkyPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=11291757 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `StaticTile`
--

DROP TABLE IF EXISTS `StaticTile`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `StaticTile` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Detector_id` int(11) NOT NULL,
  `FieldName` varchar(45) NOT NULL,
  `RA` double NOT NULL,
  `_Dec` double NOT NULL,
  `Coord` point NOT NULL,
  `Poly` multipolygon NOT NULL,
  `EBV` double NOT NULL,
  `N128_SkyPixel_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_StaticTile_Detector_Detector_id_idx` (`Detector_id`),
  KEY `FK_StaticTile_SkyPixel_N128_SkyPixel_id_idx` (`N128_SkyPixel_id`),
  CONSTRAINT `FK_StaticTile_Detector_Detector_id` FOREIGN KEY (`Detector_id`) REFERENCES `Detector` (`id`),
  CONSTRAINT `FK_StaticTile_SkyPixel_N128_SkyPixel_id` FOREIGN KEY (`N128_SkyPixel_id`) REFERENCES `SkyPixel` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=512553 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `StaticTile_HealpixPixel`
--

DROP TABLE IF EXISTS `StaticTile_HealpixPixel`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
 SET character_set_client = utf8mb4 ;
CREATE TABLE `StaticTile_HealpixPixel` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `StaticTile_id` int(11) NOT NULL,
  `HealpixPixel_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_StaticTile_HealpixPixel_StaticTile_id_idx` (`StaticTile_id`),
  KEY `FK_StaticTile_HealpixPixel_HealpixPixel_id_idx` (`HealpixPixel_id`),
  CONSTRAINT `FK_StaticTile_HealpixPixel_HealpixPixel_id` FOREIGN KEY (`HealpixPixel_id`) REFERENCES `HealpixPixel` (`id`),
  CONSTRAINT `FK_StaticTile_HealpixPixel_StaticTile_id` FOREIGN KEY (`StaticTile_id`) REFERENCES `StaticTile` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=125196860 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2022-08-16 15:42:49
